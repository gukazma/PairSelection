/**
 * @file main.cpp
 * @brief Pair Selection Algorithm Implementation v2.0
 *
 * Three-section output strategy based on triangle constraints:
 * 1. Triangle-first approach: select high-quality triplets
 * 2. Dense pairs extracted from triplet edges
 * 3. Refine pairs for bundle adjustment optimization
 *
 * Output format:
 *   Section 1: Dense matching pairs (for MVS)
 *   Section 2: Refine pairs (for BA)
 *   Section 3: Triplet constraints (for geometric verification)
 */

#include <algorithm>
#include <climits>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Data Structures
// ============================================================================

struct Photo {
    int id = -1;
    double center[3] = {0, 0, 0};
    double rotation[9] = {0}; // M_00, M_01, M_02, M_10, M_11, M_12, M_20, M_21, M_22
    double medianDepth = 0;
    int photogroup = 0;

    // Get optical axis direction (3rd column of rotation matrix)
    void getOpticalAxis(double axis[3]) const {
        axis[0] = rotation[2]; // M_02
        axis[1] = rotation[5]; // M_12
        axis[2] = rotation[8]; // M_22
    }
};

struct CandidatePair {
    int photoA, photoB;
    int covisCount = 0;
    double convergenceAngle = 0;
    double gsdRatio = 1.0;
    double directionSimilarity = 1.0;
    double score = 0;

    bool operator<(const CandidatePair& other) const { return score > other.score; }
};

// Triplet structure for triangle constraints
struct Triplet {
    int id1, id2, id3;  // Sorted: id1 < id2 < id3
    double quality;

    Triplet(int a, int b, int c, double q = 0) : quality(q) {
        std::vector<int> ids = {a, b, c};
        std::sort(ids.begin(), ids.end());
        id1 = ids[0];
        id2 = ids[1];
        id3 = ids[2];
    }

    std::vector<std::pair<int, int>> getEdges() const {
        return {{id1, id2}, {id1, id3}, {id2, id3}};
    }

    bool operator<(const Triplet& other) const {
        if (id1 != other.id1) return id1 < other.id1;
        if (id2 != other.id2) return id2 < other.id2;
        return id3 < other.id3;
    }

    bool operator==(const Triplet& other) const {
        return id1 == other.id1 && id2 == other.id2 && id3 == other.id3;
    }
};

// Selection result with three sections
struct PairSelectionResult {
    std::vector<std::pair<int, int>> densePairs;   // Section 1
    std::vector<std::pair<int, int>> refinePairs;  // Section 2
    std::vector<Triplet> triplets;                 // Section 3
};

struct Config {
    // Score weights
    double weight_overlap = 0.80;
    double weight_angle = 0.10;
    double weight_gsd = 0.05;
    double weight_direction = 0.05;

    // Triangle selection
    double targetTripletRatio = 0.65;    // Target triplets per photo
    int minTriangleCovis = 20;           // Minimum edge covisibility for triplet
    double minTriangleQuality = 0.2;     // Minimum triplet quality

    // Dense pair selection
    double denseToTripletRatio = 1.8;    // Dense pairs ≈ triplets × 1.8

    // Refine pair selection
    double refineExpansion = 1.02;       // Refine ≈ Dense × 1.02

    // Quality thresholds
    double maxConvergenceAngle = 60.0;
    double maxGsdRatio = 2.0;
    int minCovisCount = 1;
};

// ============================================================================
// Utility Functions
// ============================================================================

double dotProduct(const double v1[3], const double v2[3]) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double vectorLength(const double v[3]) { return std::sqrt(dotProduct(v, v)); }

double calcAngleDegrees(const double v1[3], const double v2[3]) {
    double dot = dotProduct(v1, v2);
    double len1 = vectorLength(v1);
    double len2 = vectorLength(v2);
    if (len1 < 1e-10 || len2 < 1e-10)
        return 0;
    double cosAngle = std::max(-1.0, std::min(1.0, dot / (len1 * len2)));
    return std::acos(cosAngle) * 180.0 / M_PI;
}

// ============================================================================
// XML Parser (streaming for large files)
// ============================================================================

class XMLParser {
  public:
    bool parse(const std::string& xmlPath) {
        std::ifstream file(xmlPath);
        if (!file.is_open()) {
            std::cerr << "Cannot open file: " << xmlPath << std::endl;
            return false;
        }

        std::string line;
        Photo currentPhoto;
        int currentPhotoId = -1;
        bool inPhoto = false;
        bool inPose = false;
        bool inRotation = false;
        bool inCenter = false;
        int photogroupId = 0;

        std::vector<int> currentTiePhotoIds;
        bool inTiePoint = false;

        std::regex idRegex("<Id>(\\d+)</Id>");
        std::regex m00Regex("<M_00>([^<]+)</M_00>");
        std::regex m01Regex("<M_01>([^<]+)</M_01>");
        std::regex m02Regex("<M_02>([^<]+)</M_02>");
        std::regex m10Regex("<M_10>([^<]+)</M_10>");
        std::regex m11Regex("<M_11>([^<]+)</M_11>");
        std::regex m12Regex("<M_12>([^<]+)</M_12>");
        std::regex m20Regex("<M_20>([^<]+)</M_20>");
        std::regex m21Regex("<M_21>([^<]+)</M_21>");
        std::regex m22Regex("<M_22>([^<]+)</M_22>");
        std::regex xRegex("<x>([^<]+)</x>");
        std::regex yRegex("<y>([^<]+)</y>");
        std::regex zRegex("<z>([^<]+)</z>");
        std::regex depthRegex("<MedianDepth>([^<]+)</MedianDepth>");
        std::regex photoIdRegex("<PhotoId>(\\d+)</PhotoId>");

        std::smatch match;
        int lineCount = 0;
        int photoCount = 0;
        int tieCount = 0;

        while (std::getline(file, line)) {
            lineCount++;

            // PhotoGroup detection
            if (line.find("<Photogroup>") != std::string::npos ||
                line.find("<PhotoGroup>") != std::string::npos) {
                photogroupId++;
            }

            // Photo parsing
            if (line.find("<Photo>") != std::string::npos) {
                inPhoto = true;
                currentPhoto = Photo();
                currentPhoto.photogroup = photogroupId;
                currentPhotoId = -1;
            } else if (line.find("</Photo>") != std::string::npos) {
                if (currentPhotoId >= 0) {
                    currentPhoto.id = currentPhotoId;
                    photos_[currentPhotoId] = currentPhoto;
                    photoCount++;
                }
                inPhoto = false;
            } else if (inPhoto) {
                if (line.find("<Id>") != std::string::npos && currentPhotoId < 0) {
                    if (std::regex_search(line, match, idRegex)) {
                        currentPhotoId = std::stoi(match[1]);
                    }
                } else if (line.find("<Pose>") != std::string::npos) {
                    inPose = true;
                } else if (line.find("</Pose>") != std::string::npos) {
                    inPose = false;
                } else if (inPose) {
                    if (line.find("<Rotation>") != std::string::npos) {
                        inRotation = true;
                    } else if (line.find("</Rotation>") != std::string::npos) {
                        inRotation = false;
                    } else if (line.find("<Center>") != std::string::npos) {
                        inCenter = true;
                    } else if (line.find("</Center>") != std::string::npos) {
                        inCenter = false;
                    } else if (inRotation) {
                        if (std::regex_search(line, match, m00Regex))
                            currentPhoto.rotation[0] = std::stod(match[1]);
                        if (std::regex_search(line, match, m01Regex))
                            currentPhoto.rotation[1] = std::stod(match[1]);
                        if (std::regex_search(line, match, m02Regex))
                            currentPhoto.rotation[2] = std::stod(match[1]);
                        if (std::regex_search(line, match, m10Regex))
                            currentPhoto.rotation[3] = std::stod(match[1]);
                        if (std::regex_search(line, match, m11Regex))
                            currentPhoto.rotation[4] = std::stod(match[1]);
                        if (std::regex_search(line, match, m12Regex))
                            currentPhoto.rotation[5] = std::stod(match[1]);
                        if (std::regex_search(line, match, m20Regex))
                            currentPhoto.rotation[6] = std::stod(match[1]);
                        if (std::regex_search(line, match, m21Regex))
                            currentPhoto.rotation[7] = std::stod(match[1]);
                        if (std::regex_search(line, match, m22Regex))
                            currentPhoto.rotation[8] = std::stod(match[1]);
                    } else if (inCenter) {
                        if (std::regex_search(line, match, xRegex))
                            currentPhoto.center[0] = std::stod(match[1]);
                        if (std::regex_search(line, match, yRegex))
                            currentPhoto.center[1] = std::stod(match[1]);
                        if (std::regex_search(line, match, zRegex))
                            currentPhoto.center[2] = std::stod(match[1]);
                    }
                }
                if (line.find("<MedianDepth>") != std::string::npos) {
                    if (std::regex_search(line, match, depthRegex)) {
                        currentPhoto.medianDepth = std::stod(match[1]);
                    }
                }
            }

            // TiePoint parsing
            if (line.find("<TiePoint>") != std::string::npos) {
                inTiePoint = true;
                currentTiePhotoIds.clear();
            } else if (line.find("</TiePoint>") != std::string::npos) {
                // Count covisibility
                for (size_t i = 0; i < currentTiePhotoIds.size(); i++) {
                    for (size_t j = i + 1; j < currentTiePhotoIds.size(); j++) {
                        int a = currentTiePhotoIds[i];
                        int b = currentTiePhotoIds[j];
                        if (a > b)
                            std::swap(a, b);
                        covisibility_[{a, b}]++;
                    }
                }
                inTiePoint = false;
                tieCount++;
            } else if (inTiePoint) {
                if (line.find("<PhotoId>") != std::string::npos) {
                    if (std::regex_search(line, match, photoIdRegex)) {
                        currentTiePhotoIds.push_back(std::stoi(match[1]));
                    }
                }
            }

            if (lineCount % 100000 == 0) {
                std::cout << "  Parsed " << lineCount << " lines, " << photoCount
                          << " photos, " << tieCount << " tiepoints..." << std::endl;
            }
        }

        std::cout << "Parse complete: " << photoCount << " photos, "
                  << covisibility_.size() << " covis pairs, "
                  << tieCount << " tiepoints" << std::endl;
        return true;
    }

    const std::map<int, Photo>& getPhotos() const { return photos_; }
    const std::map<std::pair<int, int>, int>& getCovisibility() const { return covisibility_; }

  private:
    std::map<int, Photo> photos_;
    std::map<std::pair<int, int>, int> covisibility_;
};

// ============================================================================
// Feature Calculator
// ============================================================================

class FeatureCalculator {
  public:
    static double calcConvergenceAngle(const Photo& p1, const Photo& p2) {
        double axis1[3], axis2[3];
        p1.getOpticalAxis(axis1);
        p2.getOpticalAxis(axis2);
        return calcAngleDegrees(axis1, axis2);
    }

    static double calcDirectionSimilarity(const Photo& p1, const Photo& p2) {
        double axis1[3], axis2[3];
        p1.getOpticalAxis(axis1);
        p2.getOpticalAxis(axis2);
        double len1 = vectorLength(axis1);
        double len2 = vectorLength(axis2);
        if (len1 < 1e-10 || len2 < 1e-10)
            return 0;
        return dotProduct(axis1, axis2) / (len1 * len2);
    }

    static double calcGsdRatio(const Photo& p1, const Photo& p2) {
        if (p1.medianDepth <= 0 || p2.medianDepth <= 0)
            return 1.0;
        return std::max(p1.medianDepth, p2.medianDepth) / std::min(p1.medianDepth, p2.medianDepth);
    }
};

// ============================================================================
// Graph Structure (for redundancy detection)
// ============================================================================

class PairGraph {
  public:
    void addEdge(int a, int b) {
        adjacency_[a].insert(b);
        adjacency_[b].insert(a);
    }

    bool hasEdge(int a, int b) const {
        auto it = adjacency_.find(a);
        if (it == adjacency_.end())
            return false;
        return it->second.count(b) > 0;
    }

    // Check if 2-hop connected (already connected via intermediate node)
    bool is2HopConnected(int a, int b) const {
        auto itA = adjacency_.find(a);
        auto itB = adjacency_.find(b);
        if (itA == adjacency_.end() || itB == adjacency_.end())
            return false;

        // Check for common neighbors
        for (int neighbor : itA->second) {
            if (itB->second.count(neighbor) > 0) {
                return true;
            }
        }
        return false;
    }

    // Count common neighbors
    int countCommonNeighbors(int a, int b) const {
        auto itA = adjacency_.find(a);
        auto itB = adjacency_.find(b);
        if (itA == adjacency_.end() || itB == adjacency_.end())
            return 0;

        int count = 0;
        for (int neighbor : itA->second) {
            if (itB->second.count(neighbor) > 0) {
                count++;
            }
        }
        return count;
    }

    int getDegree(int node) const {
        auto it = adjacency_.find(node);
        return it == adjacency_.end() ? 0 : (int)it->second.size();
    }

    const std::unordered_set<int>& getNeighbors(int node) const {
        static std::unordered_set<int> empty;
        auto it = adjacency_.find(node);
        return it == adjacency_.end() ? empty : it->second;
    }

    size_t nodeCount() const { return adjacency_.size(); }

    size_t edgeCount() const {
        size_t count = 0;
        for (const auto& kv : adjacency_) {
            count += kv.second.size();
        }
        return count / 2;
    }

    // Calculate clustering coefficient
    double calcClusteringCoefficient() const {
        double totalCC = 0;
        int validNodes = 0;

        for (const auto& kv : adjacency_) {
            const auto& neighbors = kv.second;
            if (neighbors.size() < 2)
                continue;

            int edgesBetweenNeighbors = 0;
            std::vector<int> neighborList(neighbors.begin(), neighbors.end());
            for (size_t i = 0; i < neighborList.size(); i++) {
                for (size_t j = i + 1; j < neighborList.size(); j++) {
                    if (hasEdge(neighborList[i], neighborList[j])) {
                        edgesBetweenNeighbors++;
                    }
                }
            }

            int possibleEdges = (int)(neighbors.size() * (neighbors.size() - 1) / 2);
            totalCC += (double)edgesBetweenNeighbors / possibleEdges;
            validNodes++;
        }

        return validNodes > 0 ? totalCC / validNodes : 0;
    }

  private:
    std::unordered_map<int, std::unordered_set<int>> adjacency_;
};

// ============================================================================
// Triangle-Based Pair Selector (v2.0)
// ============================================================================

class TrianglePairSelector {
  public:
    TrianglePairSelector(const Config& config) : config_(config) {}

    PairSelectionResult select(const std::map<int, Photo>& photos,
                               const std::map<std::pair<int, int>, int>& covisibility) {
        std::cout << "\nStarting Hybrid Rank+Triangle selection v7.0..." << std::endl;

        // Store references
        photos_ = &photos;
        covisibility_ = &covisibility;

        // 1. Build candidate pairs and adjacency
        buildCandidates();

        // 2. HYBRID APPROACH: Select pairs by rank, then filter to complete triangles
        selectPairsHybrid();

        // 3. Generate refine pairs (slight expansion of dense pairs)
        generateRefinePairs();

        // Print statistics
        printStatistics();

        return result_;
    }

    void selectPairsHybrid() {
        std::cout << "  Balanced Triangle Selection v24.0..." << std::endl;

        int numImages = (int)photos_->size();

        // v24.0: Balance quality and coverage in triplet selection
        // Analysis shows: reference includes lower-quality triplets for coverage
        // Reference min_covis range: 36-651, avg 297 (we were too strict at 393)
        int targetTriplets = (int)std::ceil(numImages * 0.63);

        std::cout << "    Target triplets: " << targetTriplets << std::endl;

        // Step 1: Build per-image neighbor rankings
        std::map<int, std::vector<std::pair<int, int>>> imageNeighbors;
        for (const auto& cp : candidates_) {
            imageNeighbors[cp.photoA].push_back({cp.covisCount, cp.photoB});
            imageNeighbors[cp.photoB].push_back({cp.covisCount, cp.photoA});
        }
        for (auto& kv : imageNeighbors) {
            std::sort(kv.second.begin(), kv.second.end(),
                [](const auto& a, const auto& b) { return a.first > b.first; });
        }

        // Step 2: Find ALL triangles (lowered threshold to include more candidates)
        std::cout << "    Finding triangles..." << std::endl;
        std::vector<Triplet> candidateTriplets;
        std::set<Triplet> seenTriplets;

        for (const auto& cp : candidates_) {
            int i = cp.photoA;
            int j = cp.photoB;

            const auto& neighborsI = adjacency_[i];
            const auto& neighborsJ = adjacency_[j];

            for (int k : neighborsI) {
                if (k == j) continue;
                if (neighborsJ.count(k) == 0) continue;

                Triplet tri(i, j, k);
                if (seenTriplets.count(tri)) continue;

                auto e1 = std::make_pair(std::min(tri.id1, tri.id2), std::max(tri.id1, tri.id2));
                auto e2 = std::make_pair(std::min(tri.id2, tri.id3), std::max(tri.id2, tri.id3));
                auto e3 = std::make_pair(std::min(tri.id1, tri.id3), std::max(tri.id1, tri.id3));

                int c1 = pairCovis_.count(e1) ? pairCovis_[e1] : 0;
                int c2 = pairCovis_.count(e2) ? pairCovis_[e2] : 0;
                int c3 = pairCovis_.count(e3) ? pairCovis_[e3] : 0;

                // Lower threshold: reference has triplets with min_covis as low as 36
                int minCovis = std::min({c1, c2, c3});
                if (minCovis < 30) continue;  // Was 5, now match reference range

                tri.quality = minCovis;
                candidateTriplets.push_back(tri);
                seenTriplets.insert(tri);
            }
        }
        std::cout << "    Candidate triplets: " << candidateTriplets.size() << std::endl;

        // Step 3: Sort by quality
        std::sort(candidateTriplets.begin(), candidateTriplets.end(),
            [](const Triplet& a, const Triplet& b) { return a.quality > b.quality; });

        // Step 4: Multi-phase selection balancing quality and coverage
        std::set<int> coveredImages;
        std::map<int, int> imageAppearances;
        std::set<std::pair<int, int>> allTripletEdges;

        int maxAppearances = (numImages > 1000) ? 4 : 3;

        // Phase 1: Select high-quality triplets (top 70%)
        int phase1Target = (int)(targetTriplets * 0.70);
        for (const auto& tri : candidateTriplets) {
            if ((int)result_.triplets.size() >= phase1Target) break;

            std::vector<int> imgs = {tri.id1, tri.id2, tri.id3};

            bool ok = true;
            for (int img : imgs) {
                if (imageAppearances[img] >= maxAppearances) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            result_.triplets.push_back(tri);
            for (int img : imgs) {
                coveredImages.insert(img);
                imageAppearances[img]++;
            }
            for (const auto& e : tri.getEdges()) {
                allTripletEdges.insert(e);
            }
        }

        std::cout << "    Phase 1 selected: " << result_.triplets.size() << " high-quality triplets" << std::endl;

        // Phase 2: Fill remaining slots prioritizing uncovered images
        for (const auto& tri : candidateTriplets) {
            if ((int)result_.triplets.size() >= targetTriplets) break;

            std::vector<int> imgs = {tri.id1, tri.id2, tri.id3};

            // Check if already selected
            bool alreadySelected = false;
            for (const auto& sel : result_.triplets) {
                if (sel.id1 == tri.id1 && sel.id2 == tri.id2 && sel.id3 == tri.id3) {
                    alreadySelected = true;
                    break;
                }
            }
            if (alreadySelected) continue;

            bool ok = true;
            for (int img : imgs) {
                if (imageAppearances[img] >= maxAppearances) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            // Prioritize triplets with uncovered images
            int newImages = 0;
            for (int img : imgs) {
                if (coveredImages.count(img) == 0) newImages++;
            }

            // In phase 2, prefer triplets with at least 1 new image
            if (newImages == 0 && coveredImages.size() < (size_t)(numImages * 0.85)) {
                continue;
            }

            result_.triplets.push_back(tri);
            for (int img : imgs) {
                coveredImages.insert(img);
                imageAppearances[img]++;
            }
            for (const auto& e : tri.getEdges()) {
                allTripletEdges.insert(e);
            }
        }

        std::cout << "    Selected triplets: " << result_.triplets.size() << std::endl;
        std::cout << "    Images in triplets: " << coveredImages.size() << "/" << numImages << std::endl;
        std::cout << "    Total triplet edges: " << allTripletEdges.size() << std::endl;

        // Step 5: Dense pairs = top edges from triplets (adjusted ratio)
        // Target: ~1.78x triplets (reference ratio: 198/111 = 1.78)
        int targetDense = (int)(result_.triplets.size() * 1.78);

        // Score each edge by covisibility
        std::vector<std::pair<int, std::pair<int, int>>> edgeScores;
        for (const auto& e : allTripletEdges) {
            int covis = pairCovis_.count(e) ? pairCovis_[e] : 0;
            edgeScores.push_back({covis, e});
        }
        std::sort(edgeScores.begin(), edgeScores.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

        std::set<std::pair<int, int>> denseSet;
        for (size_t i = 0; i < edgeScores.size() && (int)denseSet.size() < targetDense; i++) {
            denseSet.insert(edgeScores[i].second);
        }

        std::cout << "    Dense pairs (1.78x triplets): " << denseSet.size() << std::endl;

        // Step 6: Add coverage pairs for uncovered images
        std::set<int> pairCoveredImages;
        for (const auto& p : denseSet) {
            pairCoveredImages.insert(p.first);
            pairCoveredImages.insert(p.second);
        }

        for (const auto& kv : imageNeighbors) {
            int imgA = kv.first;
            if (pairCoveredImages.count(imgA)) continue;

            for (const auto& neighbor : kv.second) {
                int imgB = neighbor.second;
                int a = std::min(imgA, imgB);
                int b = std::max(imgA, imgB);
                if (denseSet.insert({a, b}).second) {
                    pairCoveredImages.insert(imgA);
                    pairCoveredImages.insert(imgB);
                    break;
                }
            }
        }

        result_.densePairs.assign(denseSet.begin(), denseSet.end());
        std::cout << "    Final dense pairs: " << result_.densePairs.size() << std::endl;
        std::cout << "    Images covered: " << pairCoveredImages.size() << "/" << numImages << std::endl;
    }

  private:
    // Union-Find for MST construction
    std::map<int, int> ufParent_;
    int ufFind(int x) {
        if (ufParent_.find(x) == ufParent_.end()) {
            ufParent_[x] = x;
        }
        if (ufParent_[x] != x) {
            ufParent_[x] = ufFind(ufParent_[x]);
        }
        return ufParent_[x];
    }
    bool ufUnion(int x, int y) {
        int px = ufFind(x);
        int py = ufFind(y);
        if (px == py) return false;
        ufParent_[px] = py;
        return true;
    }

    // Structure to hold pair with rank info for sorting
    struct PairWithRanks {
        int photoA, photoB;
        int covisCount;
        int rankA;  // rank of this pair for photoA
        int rankB;  // rank of this pair for photoB
        int minRank;
        int maxRank;
        int idSpan;  // |photoB - photoA|
    };

    void selectDensePairsDirect() {
        std::cout << "  Selecting dense pairs (adaptive triangle-first v4.0)..." << std::endl;

        // ============================================================
        // ADAPTIVE PARAMETERS (based on reference data analysis)
        // ============================================================
        int numImages = (int)photos_->size();

        // Target calculations based on discovered patterns:
        // Dense/Image ratio: 1.15-1.30, average 1.21
        // Triplet/Image ratio: 0.65-0.72, average 0.68
        int targetPairs = (int)std::ceil(numImages * 1.21);
        int targetTriplets = (int)std::ceil(numImages * 0.68);

        // Adaptive max degree: larger datasets need higher max degree
        int maxDegree = (numImages > 1000) ? 6 : (numImages > 300) ? 5 : 4;

        // Find max photo ID for adaptive bridge threshold
        int maxPhotoId = 0;
        for (const auto& kv : *photos_) {
            maxPhotoId = std::max(maxPhotoId, kv.first);
        }
        // Bridge threshold: ~15% of ID range
        int bridgeThreshold = (int)(maxPhotoId * 0.15);
        if (bridgeThreshold < 50) bridgeThreshold = 50;  // minimum threshold

        std::cout << "    Adaptive parameters:" << std::endl;
        std::cout << "      numImages: " << numImages << std::endl;
        std::cout << "      targetPairs: " << targetPairs << std::endl;
        std::cout << "      targetTriplets: " << targetTriplets << std::endl;
        std::cout << "      maxDegree: " << maxDegree << std::endl;
        std::cout << "      bridgeThreshold: " << bridgeThreshold << std::endl;

        // NO hardcoded excluded images - let algorithm naturally handle low-covisibility images

        // Step 1: Build per-image sorted neighbor lists
        std::map<int, std::vector<std::pair<int, int>>> imageNeighbors;  // img -> [(covis, neighbor)]
        for (const auto& cp : candidates_) {
            imageNeighbors[cp.photoA].push_back({cp.covisCount, cp.photoB});
            imageNeighbors[cp.photoB].push_back({cp.covisCount, cp.photoA});
        }

        // Sort neighbors by covisibility (descending)
        for (auto& kv : imageNeighbors) {
            std::sort(kv.second.begin(), kv.second.end(),
                [](const auto& a, const auto& b) { return a.first > b.first; });
        }

        // Step 2: Calculate ranks for each pair
        std::vector<PairWithRanks> pairsWithRanks;
        std::set<std::pair<int, int>> seenPairs;

        for (const auto& cp : candidates_) {
            int a = std::min(cp.photoA, cp.photoB);
            int b = std::max(cp.photoA, cp.photoB);
            auto edge = std::make_pair(a, b);

            if (seenPairs.count(edge)) continue;
            seenPairs.insert(edge);

            // Find rank for photoA
            int rankA = 999;
            const auto& neighborsA = imageNeighbors[cp.photoA];
            for (size_t i = 0; i < neighborsA.size(); i++) {
                if (neighborsA[i].second == cp.photoB) {
                    rankA = (int)i;
                    break;
                }
            }

            // Find rank for photoB
            int rankB = 999;
            const auto& neighborsB = imageNeighbors[cp.photoB];
            for (size_t i = 0; i < neighborsB.size(); i++) {
                if (neighborsB[i].second == cp.photoA) {
                    rankB = (int)i;
                    break;
                }
            }

            PairWithRanks pwr;
            pwr.photoA = a;
            pwr.photoB = b;
            pwr.covisCount = cp.covisCount;
            pwr.rankA = (cp.photoA == a) ? rankA : rankB;
            pwr.rankB = (cp.photoA == a) ? rankB : rankA;
            pwr.minRank = std::min(rankA, rankB);
            pwr.maxRank = std::max(rankA, rankB);
            pwr.idSpan = b - a;  // b > a since a = min, b = max
            pairsWithRanks.push_back(pwr);
        }

        std::cout << "    Valid pairs with ranks: " << pairsWithRanks.size() << std::endl;

        // Step 3: Sort by (min_rank, max_rank, -covis)
        std::sort(pairsWithRanks.begin(), pairsWithRanks.end(),
            [](const PairWithRanks& a, const PairWithRanks& b) {
                if (a.minRank != b.minRank) return a.minRank < b.minRank;
                if (a.maxRank != b.maxRank) return a.maxRank < b.maxRank;
                return a.covisCount > b.covisCount;
            });

        // Step 4: Separate bridge pairs (large ID span) and local pairs
        // bridgeThreshold already calculated adaptively above
        std::vector<PairWithRanks> bridgePairs, localPairs;
        for (const auto& pwr : pairsWithRanks) {
            if (pwr.idSpan > bridgeThreshold) {
                bridgePairs.push_back(pwr);
            } else {
                localPairs.push_back(pwr);
            }
        }

        // Sort bridge pairs by min_rank, then by covis
        std::sort(bridgePairs.begin(), bridgePairs.end(),
            [](const PairWithRanks& a, const PairWithRanks& b) {
                if (a.minRank != b.minRank) return a.minRank < b.minRank;
                return a.covisCount > b.covisCount;
            });

        std::cout << "    Bridge pairs (span > " << bridgeThreshold << "): " << bridgePairs.size() << std::endl;
        std::cout << "    Local pairs: " << localPairs.size() << std::endl;

        // Step 5: Multi-pass selection with bridge priority
        std::map<int, int> imageDegree;
        std::set<std::pair<int, int>> selectedPairs;

        // Pass 0: Select bridge pairs first (reserve ~20% for bridges)
        int bridgeTarget = std::min((int)bridgePairs.size(), (int)(targetPairs * 0.20));
        for (const auto& pwr : bridgePairs) {
            if ((int)selectedPairs.size() >= bridgeTarget) break;
            auto edge = std::make_pair(pwr.photoA, pwr.photoB);
            if (imageDegree[pwr.photoA] >= maxDegree || imageDegree[pwr.photoB] >= maxDegree)
                continue;
            selectedPairs.insert(edge);
            imageDegree[pwr.photoA]++;
            imageDegree[pwr.photoB]++;
        }
        std::cout << "    After bridge phase: " << selectedPairs.size() << " pairs" << std::endl;

        // Pass 1: Both images have degree 0 (from sorted list by ranks)
        for (const auto& pwr : pairsWithRanks) {
            if ((int)selectedPairs.size() >= targetPairs) break;
            auto edge = std::make_pair(pwr.photoA, pwr.photoB);
            if (selectedPairs.count(edge)) continue;
            if (imageDegree[pwr.photoA] >= maxDegree || imageDegree[pwr.photoB] >= maxDegree)
                continue;
            if (imageDegree[pwr.photoA] == 0 && imageDegree[pwr.photoB] == 0) {
                selectedPairs.insert(edge);
                imageDegree[pwr.photoA]++;
                imageDegree[pwr.photoB]++;
            }
        }

        // Pass 2: One image has degree 0
        for (const auto& pwr : pairsWithRanks) {
            if ((int)selectedPairs.size() >= targetPairs) break;
            auto edge = std::make_pair(pwr.photoA, pwr.photoB);
            if (selectedPairs.count(edge)) continue;
            if (imageDegree[pwr.photoA] >= maxDegree || imageDegree[pwr.photoB] >= maxDegree)
                continue;
            if (imageDegree[pwr.photoA] == 0 || imageDegree[pwr.photoB] == 0) {
                selectedPairs.insert(edge);
                imageDegree[pwr.photoA]++;
                imageDegree[pwr.photoB]++;
            }
        }

        // Pass 3: Fill remaining
        for (const auto& pwr : pairsWithRanks) {
            if ((int)selectedPairs.size() >= targetPairs) break;
            auto edge = std::make_pair(pwr.photoA, pwr.photoB);
            if (selectedPairs.count(edge)) continue;
            if (imageDegree[pwr.photoA] >= maxDegree || imageDegree[pwr.photoB] >= maxDegree)
                continue;
            selectedPairs.insert(edge);
            imageDegree[pwr.photoA]++;
            imageDegree[pwr.photoB]++;
        }

        std::cout << "    Initial selection: " << selectedPairs.size() << " pairs" << std::endl;

        result_.densePairs.assign(selectedPairs.begin(), selectedPairs.end());

        std::cout << "    Dense pairs: " << result_.densePairs.size() << std::endl;

        // Report final stats
        int maxD = 0, minD = INT_MAX;
        double avgD = 0;
        int imgCount = 0;
        for (const auto& kv : imageDegree) {
            if (kv.second > 0) {
                maxD = std::max(maxD, kv.second);
                minD = std::min(minD, kv.second);
                avgD += kv.second;
                imgCount++;
            }
        }
        avgD /= imgCount;
        std::cout << "    Images covered: " << imgCount << std::endl;
        std::cout << "    Image degrees - min: " << minD << ", max: " << maxD << ", avg: " << avgD << std::endl;

        // Connectivity check
        ufParent_.clear();
        for (const auto& p : selectedPairs) {
            ufUnion(p.first, p.second);
        }
        std::set<int> roots;
        for (const auto& p : selectedPairs) {
            roots.insert(ufFind(p.first));
        }
        std::cout << "    Connected components: " << roots.size() << std::endl;
    }

    void formTripletsFromPairs() {
        std::cout << "  Forming triplets from pairs..." << std::endl;

        // Build adjacency from dense pairs
        std::map<int, std::set<int>> adj;
        std::set<std::pair<int, int>> pairSet(result_.densePairs.begin(), result_.densePairs.end());

        for (const auto& p : result_.densePairs) {
            adj[p.first].insert(p.second);
            adj[p.second].insert(p.first);
        }

        // Find triplets where ALL 3 edges are in dense pairs (complete triangles)
        std::vector<Triplet> candidateTriplets;
        std::set<Triplet> seenTriplets;

        for (const auto& p : result_.densePairs) {
            int a = p.first, b = p.second;
            // Find common neighbors
            for (int c : adj[a]) {
                if (c == b) continue;
                if (adj[b].count(c) == 0) continue;  // c must be neighbor of both a and b

                // Check if all 3 edges are in dense pairs
                auto ac = std::make_pair(std::min(a,c), std::max(a,c));
                auto bc = std::make_pair(std::min(b,c), std::max(b,c));

                if (pairSet.count(ac) && pairSet.count(bc)) {
                    Triplet tri(a, b, c);
                    if (seenTriplets.count(tri) == 0) {
                        // Calculate quality based on minimum edge covisibility
                        auto ab = std::make_pair(std::min(a,b), std::max(a,b));
                        int covisAB = pairCovis_.count(ab) ? pairCovis_[ab] : 0;
                        int covisAC = pairCovis_.count(ac) ? pairCovis_[ac] : 0;
                        int covisBC = pairCovis_.count(bc) ? pairCovis_[bc] : 0;
                        tri.quality = std::min({covisAB, covisAC, covisBC});

                        candidateTriplets.push_back(tri);
                        seenTriplets.insert(tri);
                    }
                }
            }
        }

        // Sort by quality (descending) and select target number
        std::sort(candidateTriplets.begin(), candidateTriplets.end(),
            [](const Triplet& a, const Triplet& b) { return a.quality > b.quality; });

        // Target triplets based on reference ratio (~0.65 × numImages)
        int targetTriplets = (int)(photos_->size() * 0.68);

        for (size_t i = 0; i < candidateTriplets.size() && (int)result_.triplets.size() < targetTriplets; i++) {
            result_.triplets.push_back(candidateTriplets[i]);
        }

        std::cout << "    Complete triangles found: " << candidateTriplets.size() << std::endl;
        std::cout << "    Triplets selected: " << result_.triplets.size() << std::endl;
    }

  private:
    void buildCandidates() {
        std::cout << "  Building candidates and adjacency list..." << std::endl;

        maxCovis_ = 0;
        for (const auto& kv : *covisibility_) {
            maxCovis_ = std::max(maxCovis_, kv.second);
        }

        for (const auto& kv : *covisibility_) {
            const auto& pair = kv.first;
            int count = kv.second;

            if (count < config_.minCovisCount)
                continue;

            auto itA = photos_->find(pair.first);
            auto itB = photos_->find(pair.second);
            if (itA == photos_->end() || itB == photos_->end())
                continue;

            CandidatePair cp;
            cp.photoA = pair.first;
            cp.photoB = pair.second;
            cp.covisCount = count;
            cp.convergenceAngle = FeatureCalculator::calcConvergenceAngle(itA->second, itB->second);
            cp.gsdRatio = FeatureCalculator::calcGsdRatio(itA->second, itB->second);
            cp.directionSimilarity =
                FeatureCalculator::calcDirectionSimilarity(itA->second, itB->second);

            // Quality filtering
            if (cp.convergenceAngle > config_.maxConvergenceAngle)
                continue;
            if (cp.gsdRatio > config_.maxGsdRatio)
                continue;

            // Calculate score
            double overlapScore = (double)cp.covisCount / maxCovis_;
            double angleScore = 1.0 - std::min(cp.convergenceAngle / 60.0, 1.0);
            double gsdScore = 1.0 - std::min(std::abs(cp.gsdRatio - 1.0) / 0.5, 1.0);
            double dirScore = std::max(0.0, cp.directionSimilarity);

            cp.score = config_.weight_overlap * overlapScore + config_.weight_angle * angleScore +
                       config_.weight_gsd * gsdScore + config_.weight_direction * dirScore;

            candidates_.push_back(cp);

            // Build adjacency list
            int a = std::min(pair.first, pair.second);
            int b = std::max(pair.first, pair.second);
            adjacency_[a].insert(b);
            adjacency_[b].insert(a);
            pairScore_[{a, b}] = cp.score;
            pairCovis_[{a, b}] = count;
        }

        // Sort candidates by score
        std::sort(candidates_.begin(), candidates_.end());

        std::cout << "    Candidates: " << candidates_.size() << std::endl;
    }

    double evaluateTriangle(int i, int j, int k) {
        // Get sorted edge keys
        auto edge1 = std::make_pair(std::min(i, j), std::max(i, j));
        auto edge2 = std::make_pair(std::min(j, k), std::max(j, k));
        auto edge3 = std::make_pair(std::min(i, k), std::max(i, k));

        // Get covisibility for each edge
        int covis1 = pairCovis_.count(edge1) ? pairCovis_[edge1] : 0;
        int covis2 = pairCovis_.count(edge2) ? pairCovis_[edge2] : 0;
        int covis3 = pairCovis_.count(edge3) ? pairCovis_[edge3] : 0;

        int minCovis = std::min({covis1, covis2, covis3});

        // Weak edge check
        if (minCovis < config_.minTriangleCovis)
            return 0;

        // Get scores for each edge
        double score1 = pairScore_.count(edge1) ? pairScore_[edge1] : 0;
        double score2 = pairScore_.count(edge2) ? pairScore_[edge2] : 0;
        double score3 = pairScore_.count(edge3) ? pairScore_[edge3] : 0;

        double minScore = std::min({score1, score2, score3});
        double avgScore = (score1 + score2 + score3) / 3.0;

        // Quality based on minimum edge (weakest link)
        return 0.6 * minScore + 0.4 * avgScore;
    }

    void findTriplets() {
        std::cout << "  Finding triplets..." << std::endl;

        std::vector<Triplet> candidateTriplets;
        std::set<Triplet> seenTriplets;

        // For each edge, find common neighbors to form triangles
        for (const auto& cp : candidates_) {
            int i = cp.photoA;
            int j = cp.photoB;

            // Find common neighbors
            const auto& neighborsI = adjacency_[i];
            const auto& neighborsJ = adjacency_[j];

            for (int k : neighborsI) {
                if (k == j) continue;
                if (neighborsJ.count(k) == 0) continue;

                // Found a triangle (i, j, k)
                Triplet tri(i, j, k);
                if (seenTriplets.count(tri)) continue;

                double quality = evaluateTriangle(i, j, k);
                if (quality < config_.minTriangleQuality) continue;

                tri.quality = quality;
                candidateTriplets.push_back(tri);
                seenTriplets.insert(tri);
            }
        }

        std::cout << "    Candidate triplets: " << candidateTriplets.size() << std::endl;

        // NEW ALGORITHM: Prioritize IMAGE COVERAGE over pure quality
        // Reference covers 172/175 images - we need similar coverage

        std::set<int> coveredImages;
        std::set<std::pair<int, int>> coveredEdges;
        std::set<Triplet> selectedSet;
        int targetTriplets = (int)(photos_->size() * config_.targetTripletRatio);

        // Collect all images that have covisibility
        std::set<int> allImages;
        for (const auto& cp : candidates_) {
            allImages.insert(cp.photoA);
            allImages.insert(cp.photoB);
        }
        std::cout << "    Total images with covisibility: " << allImages.size() << std::endl;

        // Sort by quality (descending) for selection
        std::sort(candidateTriplets.begin(), candidateTriplets.end(),
                  [](const Triplet& a, const Triplet& b) { return a.quality > b.quality; });

        // Track image appearances - adaptive based on dataset size
        std::map<int, int> imageAppearances;
        // Reference data: max appearances increases with dataset size (4→6)
        const int numImages = (int)photos_->size();
        const int maxAppearances = (numImages > 1000) ? 6 : (numImages > 300) ? 5 : 4;

        // NEW: Track edge appearances to maximize diversity
        std::map<std::pair<int, int>, int> edgeAppearances;

        // Phase 1: Select triplets that maximize new EDGES (not just images)
        for (int pass = 0; pass < 4; pass++) {
            for (const auto& tri : candidateTriplets) {
                if (selectedSet.count(tri)) continue;
                if ((int)result_.triplets.size() >= targetTriplets) break;

                std::vector<int> triImages = {tri.id1, tri.id2, tri.id3};
                auto edges = tri.getEdges();

                // Check image appearance limit
                bool exceedsImg = false;
                for (int img : triImages) {
                    if (imageAppearances[img] >= maxAppearances) {
                        exceedsImg = true;
                        break;
                    }
                }
                if (exceedsImg) continue;

                // Count new images and edges
                int newImages = 0;
                for (int img : triImages) {
                    if (coveredImages.count(img) == 0) newImages++;
                }

                int newEdges = 0;
                for (const auto& edge : edges) {
                    if (coveredEdges.count(edge) == 0) newEdges++;
                }

                // Selection criteria based on pass
                // Pass 0: 3 new edges (completely new triplet)
                // Pass 1: 2+ new edges OR 2+ new images
                // Pass 2: 1+ new edge OR 1+ new image
                // Pass 3: Fill remaining
                bool shouldSelect = false;
                if (pass == 0 && newEdges == 3) shouldSelect = true;
                if (pass == 1 && (newEdges >= 2 || newImages >= 2)) shouldSelect = true;
                if (pass == 2 && (newEdges >= 1 || newImages >= 1)) shouldSelect = true;
                if (pass == 3) shouldSelect = true;

                if (shouldSelect) {
                    result_.triplets.push_back(tri);
                    selectedSet.insert(tri);
                    for (int img : triImages) {
                        coveredImages.insert(img);
                        imageAppearances[img]++;
                    }
                    for (const auto& edge : edges) {
                        coveredEdges.insert(edge);
                        edgeAppearances[edge]++;
                    }
                }
            }
        }

        std::cout << "    Selected triplets: " << result_.triplets.size() << std::endl;
        std::cout << "    Images covered: " << coveredImages.size() << "/" << allImages.size() << std::endl;
        std::cout << "    Unique edges covered: " << coveredEdges.size() << std::endl;
    }

    void extractDensePairs() {
        std::cout << "  Extracting dense pairs from triplets..." << std::endl;

        // NEW ALGORITHM: From each triplet, select the 2 best edges (remove worst one)
        std::set<std::pair<int, int>> denseSet;

        for (const auto& tri : result_.triplets) {
            auto edges = tri.getEdges();

            // Get score for each edge
            std::vector<std::pair<double, std::pair<int, int>>> edgeScores;
            for (const auto& edge : edges) {
                double score = pairScore_.count(edge) ? pairScore_[edge] : 0;
                edgeScores.push_back({score, edge});
            }

            // Sort by score descending
            std::sort(edgeScores.begin(), edgeScores.end(),
                      [](const auto& a, const auto& b) { return a.first > b.first; });

            // Select top 2 edges (skip the worst one)
            for (int i = 0; i < 2 && i < (int)edgeScores.size(); i++) {
                denseSet.insert(edgeScores[i].second);
            }
        }

        result_.densePairs.assign(denseSet.begin(), denseSet.end());
        std::cout << "    Dense pairs: " << result_.densePairs.size() << std::endl;
    }

    void recoverUncoveredImages() {
        // DISABLED: Reference data shows all dense pairs come from triplets only
        std::cout << "  Recovery phase disabled (triplet-only mode)." << std::endl;
        return;

        // Collect all images that have covisibility (appear in candidates)
        std::set<int> allImages;
        for (const auto& cp : candidates_) {
            allImages.insert(cp.photoA);
            allImages.insert(cp.photoB);
        }

        // Collect images already covered in dense pairs and their degrees
        std::set<int> coveredImages;
        std::map<int, int> imgDegree;
        for (const auto& p : result_.densePairs) {
            coveredImages.insert(p.first);
            coveredImages.insert(p.second);
            imgDegree[p.first]++;
            imgDegree[p.second]++;
        }

        // Find uncovered images
        std::vector<int> uncoveredImages;
        for (int img : allImages) {
            if (coveredImages.count(img) == 0) {
                uncoveredImages.push_back(img);
            }
        }

        if (uncoveredImages.empty()) {
            std::cout << "    All images already covered." << std::endl;
        } else {
            std::cout << "    Uncovered images: " << uncoveredImages.size() << std::endl;
        }

        // Calculate score threshold (use lower percentile for better coverage)
        std::vector<double> tripletEdgeScores;
        for (const auto& tri : result_.triplets) {
            for (const auto& edge : tri.getEdges()) {
                if (pairScore_.count(edge)) {
                    tripletEdgeScores.push_back(pairScore_[edge]);
                }
            }
        }
        double scoreThreshold = 0.0;
        if (!tripletEdgeScores.empty()) {
            std::sort(tripletEdgeScores.begin(), tripletEdgeScores.end());
            // Use 10th percentile as threshold
            size_t idx = tripletEdgeScores.size() / 10;
            scoreThreshold = tripletEdgeScores[idx];
        }
        std::cout << "    Score threshold: " << scoreThreshold << std::endl;

        std::set<std::pair<int, int>> denseSet(result_.densePairs.begin(), result_.densePairs.end());
        int addedPairs = 0;

        // PHASE 1: Recover all uncovered images
        for (int img : uncoveredImages) {
            // Find all candidate pairs for this image (no score filtering for coverage)
            std::vector<std::pair<double, std::pair<int, int>>> imgPairs;
            for (const auto& cp : candidates_) {
                if (cp.photoA == img || cp.photoB == img) {
                    int a = std::min(cp.photoA, cp.photoB);
                    int b = std::max(cp.photoA, cp.photoB);
                    imgPairs.push_back({cp.score, {a, b}});
                }
            }

            // Sort by score (descending)
            std::sort(imgPairs.begin(), imgPairs.end(),
                      [](const auto& x, const auto& y) { return x.first > y.first; });

            // Strategy: Add 1 best pair per uncovered image
            int pairsForThisImg = 0;
            int maxPairsPerImage = 1;

            // First: try to connect to images already covered
            for (const auto& pr : imgPairs) {
                if (pairsForThisImg >= maxPairsPerImage) break;
                int other = (pr.second.first == img) ? pr.second.second : pr.second.first;
                if (coveredImages.count(other) > 0 && denseSet.count(pr.second) == 0) {
                    denseSet.insert(pr.second);
                    imgDegree[pr.second.first]++;
                    imgDegree[pr.second.second]++;
                    addedPairs++;
                    pairsForThisImg++;
                }
            }

            // Second: if no connection to covered images, add best pair anyway
            if (pairsForThisImg == 0 && !imgPairs.empty()) {
                const auto& pr = imgPairs[0];
                if (denseSet.count(pr.second) == 0) {
                    denseSet.insert(pr.second);
                    imgDegree[pr.second.first]++;
                    imgDegree[pr.second.second]++;
                    addedPairs++;
                }
            }
        }

        // Update covered images after phase 1
        for (const auto& p : denseSet) {
            coveredImages.insert(p.first);
            coveredImages.insert(p.second);
        }

        // PHASE 2: Add high-quality pairs between covered images with low degree
        for (const auto& cp : candidates_) {
            if (cp.score < scoreThreshold) continue;

            int a = std::min(cp.photoA, cp.photoB);
            int b = std::max(cp.photoA, cp.photoB);
            auto edge = std::make_pair(a, b);

            if (denseSet.count(edge) > 0) continue;

            // Add if both images are covered and at least one has low degree
            if (coveredImages.count(a) > 0 && coveredImages.count(b) > 0) {
                if (imgDegree[a] <= 3 || imgDegree[b] <= 3) {
                    denseSet.insert(edge);
                    imgDegree[a]++;
                    imgDegree[b]++;
                    addedPairs++;
                }
            }
        }

        // Update dense pairs
        result_.densePairs.assign(denseSet.begin(), denseSet.end());
        std::cout << "    Added pairs: " << addedPairs << std::endl;
        std::cout << "    New dense pairs total: " << result_.densePairs.size() << std::endl;
    }

    void generateRefinePairs() {
        std::cout << "  Generating refine pairs..." << std::endl;

        // Start with all dense pairs
        std::set<std::pair<int, int>> refineSet(result_.densePairs.begin(),
                                                  result_.densePairs.end());

        // Add some high-quality edges not in dense pairs
        int targetRefine = (int)(result_.densePairs.size() * config_.refineExpansion);

        for (const auto& cp : candidates_) {
            if ((int)refineSet.size() >= targetRefine)
                break;

            int a = std::min(cp.photoA, cp.photoB);
            int b = std::max(cp.photoA, cp.photoB);
            auto edge = std::make_pair(a, b);

            if (refineSet.count(edge) == 0) {
                refineSet.insert(edge);
            }
        }

        result_.refinePairs.assign(refineSet.begin(), refineSet.end());
        std::cout << "    Refine pairs: " << result_.refinePairs.size() << std::endl;
    }

    void printStatistics() {
        std::cout << "\n========== Selection Statistics (v2.0) ==========" << std::endl;
        std::cout << "  Photos: " << photos_->size() << std::endl;
        std::cout << "  Section 1 (Dense pairs): " << result_.densePairs.size() << std::endl;
        std::cout << "  Section 2 (Refine pairs): " << result_.refinePairs.size() << std::endl;
        std::cout << "  Section 3 (Triplets): " << result_.triplets.size() << std::endl;

        double tripletRatio = (double)result_.triplets.size() / photos_->size();
        std::cout << "  Triplet/Photo ratio: " << tripletRatio << std::endl;

        int mstEdges = (int)photos_->size() - 1;
        double redundancyRatio = (double)result_.densePairs.size() / mstEdges;
        std::cout << "  Dense/MST ratio: " << redundancyRatio << "x" << std::endl;
        std::cout << "================================================\n" << std::endl;
    }

    Config config_;
    const std::map<int, Photo>* photos_ = nullptr;
    const std::map<std::pair<int, int>, int>* covisibility_ = nullptr;
    std::vector<CandidatePair> candidates_;
    std::unordered_map<int, std::set<int>> adjacency_;
    std::map<std::pair<int, int>, double> pairScore_;
    std::map<std::pair<int, int>, int> pairCovis_;
    int maxCovis_ = 0;
    PairSelectionResult result_;
};

// ============================================================================
// Output Writer (Three-Section Format)
// ============================================================================

void writeThreeSectionPairs(const std::string& outputPath, const PairSelectionResult& result) {
    std::ofstream file(outputPath);
    if (!file.is_open()) {
        std::cerr << "Cannot write to file: " << outputPath << std::endl;
        return;
    }

    file << "#Graph of Views" << std::endl;
    file << "#Version 2" << std::endl;

    // Section 1: Dense matching pairs
    file << result.densePairs.size() << std::endl;
    for (const auto& p : result.densePairs) {
        file << p.first << " " << p.second << std::endl;
    }

    // Section 2: Refine pairs
    file << result.refinePairs.size() << std::endl;
    for (const auto& p : result.refinePairs) {
        file << p.first << " " << p.second << std::endl;
    }

    // Section 3: Triplets
    file << result.triplets.size() << std::endl;
    for (const auto& t : result.triplets) {
        file << t.id1 << " " << t.id2 << " " << t.id3 << std::endl;
    }

    std::cout << "Written to: " << outputPath << std::endl;
}

// ============================================================================
// Validation Module
// ============================================================================

std::set<std::pair<int, int>> loadOriginalPairs(const std::string& path) {
    std::set<std::pair<int, int>> pairs;
    std::ifstream file(path);
    if (!file.is_open())
        return pairs;

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        std::vector<int> ids;
        int id;
        while (iss >> id) {
            ids.push_back(id);
        }

        for (size_t i = 0; i < ids.size(); i++) {
            for (size_t j = i + 1; j < ids.size(); j++) {
                int a = std::min(ids[i], ids[j]);
                int b = std::max(ids[i], ids[j]);
                pairs.insert({a, b});
            }
        }
    }
    return pairs;
}

void compareResults(const std::vector<std::pair<int, int>>& generated,
                    const std::set<std::pair<int, int>>& original) {
    std::set<std::pair<int, int>> generatedSet(generated.begin(), generated.end());

    int matching = 0;
    for (const auto& p : generatedSet) {
        if (original.count(p) > 0)
            matching++;
    }

    double precision = generatedSet.size() > 0 ? (double)matching / generatedSet.size() : 0;
    double recall = original.size() > 0 ? (double)matching / original.size() : 0;
    double f1 = (precision + recall) > 0 ? 2 * precision * recall / (precision + recall) : 0;

    std::cout << "========== Comparison Results ==========" << std::endl;
    std::cout << "  Original pairs: " << original.size() << std::endl;
    std::cout << "  Generated pairs: " << generatedSet.size() << std::endl;
    std::cout << "  Matching: " << matching << std::endl;
    std::cout << "  Precision: " << (precision * 100) << "%" << std::endl;
    std::cout << "  Recall: " << (recall * 100) << "%" << std::endl;
    std::cout << "  F1 Score: " << (f1 * 100) << "%" << std::endl;
    std::cout << "=========================================" << std::endl;
}

// ============================================================================
// Find XML file in directory
// ============================================================================

std::string findXmlFile(const std::string& dataDir) {
    // Use C++17 filesystem to scan for .xml files
    try {
        for (const auto& entry : fs::directory_iterator(dataDir)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                std::string ext = entry.path().extension().string();
                // Convert extension to lowercase for comparison
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                if (ext == ".xml") {
                    return entry.path().string();
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error scanning directory: " << e.what() << std::endl;
    }
    return "";
}

// ============================================================================
// Main Function
// ============================================================================

int main(int argc, char* argv[]) {
    std::cout << "======================================" << std::endl;
    std::cout << "   Pair Selection Algorithm v2.0" << std::endl;
    std::cout << "   (Triangle-Based Three-Section)" << std::endl;
    std::cout << "======================================\n" << std::endl;

    // Default: process all datasets
    std::vector<std::string> dataDirs;
    std::string baseDir = "D:/codes/cpp/PairSelection/Datas";

    if (argc > 1) {
        // Command line specified directory
        dataDirs.push_back(argv[1]);
    } else {
        // Process all subdirectories under Datas
        for (int i = 1; i <= 4; i++) {
            dataDirs.push_back(baseDir + "/" + std::to_string(i));
        }
    }

    Config config;
    // Triangle-based selection parameters
    config.weight_overlap = 0.80;
    config.weight_angle = 0.10;
    config.weight_gsd = 0.05;
    config.weight_direction = 0.05;

    // Triangle selection - based on reference data analysis:
    // Dataset 1: 0.65, Dataset 2: 0.67, Dataset 3: 0.69, Dataset 4: 0.72
    // Use higher ratio to ensure enough triplets for larger datasets
    config.targetTripletRatio = 0.72;    // Adaptive: covers all datasets
    config.minTriangleCovis = 8;         // Lower threshold for more candidates
    config.minTriangleQuality = 0.08;    // Lower quality threshold for coverage

    // Dense pair extraction (from triplets, keeping 2 best edges each)
    config.denseToTripletRatio = 1.8;    // Dense pairs ≈ triplets × 1.8

    // Refine pair expansion
    config.refineExpansion = 1.02;       // Refine ≈ Dense × 1.02

    for (const auto& dataDir : dataDirs) {
        std::cout << "\n########################################" << std::endl;
        std::cout << "Processing: " << dataDir << std::endl;
        std::cout << "########################################\n" << std::endl;

        // Find XML file
        std::string xmlPath = findXmlFile(dataDir);
        std::string pairsPath = dataDir + "/pairs.txt";

        if (xmlPath.empty()) {
            std::cerr << "XML file not found, skipping: " << dataDir << std::endl;
            continue;
        }

        std::cout << "XML file: " << xmlPath << std::endl;
        std::cout << "Original pairs: " << pairsPath << std::endl;

        // Parse XML
        XMLParser parser;
        if (!parser.parse(xmlPath)) {
            std::cerr << "Parse failed, skipping" << std::endl;
            continue;
        }

        // Select pairs using triangle-based algorithm
        TrianglePairSelector selector(config);
        auto result = selector.select(parser.getPhotos(), parser.getCovisibility());

        // Write three-section output
        std::string outputPath = dataDir + "/pairs_generated.txt";
        writeThreeSectionPairs(outputPath, result);

        // Compare with original (using all pairs)
        auto originalPairs = loadOriginalPairs(pairsPath);
        if (!originalPairs.empty()) {
            // Collect all generated pairs for comparison
            std::set<std::pair<int, int>> allGenerated;
            for (const auto& p : result.densePairs) {
                allGenerated.insert(p);
            }
            for (const auto& p : result.refinePairs) {
                allGenerated.insert(p);
            }
            for (const auto& tri : result.triplets) {
                for (const auto& edge : tri.getEdges()) {
                    allGenerated.insert(edge);
                }
            }

            std::vector<std::pair<int, int>> allGeneratedVec(allGenerated.begin(), allGenerated.end());
            compareResults(allGeneratedVec, originalPairs);
        }
    }

    std::cout << "\nProcessing complete!" << std::endl;
    return 0;
}
