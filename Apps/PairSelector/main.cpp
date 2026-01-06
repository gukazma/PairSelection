/**
 * @file main.cpp
 * @brief Pair Selection Algorithm Implementation
 *
 * Based on reverse engineering analysis, implements pair selection strategy:
 * 1. Smart redundancy avoidance (2-hop connectivity check)
 * 2. Composite scoring (overlap + angle + GSD ratio + direction similarity)
 * 3. Triangle closure design
 * 4. Precise redundancy control (~2x MST)
 */

#include <algorithm>
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

struct Config {
    // Score weights - adjusted to match reference
    double weight_overlap = 0.80;    // Higher overlap priority
    double weight_angle = 0.05;      // Minimal angle penalty
    double weight_gsd = 0.10;
    double weight_direction = 0.05;

    // Selection parameters
    int minNeighborsPerPhoto = 2;
    int maxNeighborsPerPhoto = 3;    // Reduced to limit pair count
    double targetRedundancyRatio = 1.30;  // Match reference ~1.25-1.42x MST

    // Quality thresholds
    double maxConvergenceAngle = 60.0;
    double maxGsdRatio = 2.0;
    int minCovisCount = 1;

    // Redundancy control
    bool enableRedundancyCheck = true;
    int maxCommonNeighborsForRedundancy = 1;  // Very strict

    // Angle diversity
    bool enableAngleDiversity = false;  // Disabled
    double angleDiversityThreshold = 10.0;
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
// Pair Selector
// ============================================================================

class PairSelector {
  public:
    PairSelector(const Config& config) : config_(config) {}

    std::vector<std::pair<int, int>> select(const std::map<int, Photo>& photos,
                                            const std::map<std::pair<int, int>, int>& covisibility) {
        std::cout << "\nStarting pair selection..." << std::endl;

        // 1. Build candidate pairs and calculate features
        std::vector<CandidatePair> candidates;
        candidates.reserve(covisibility.size());

        int maxCovis = 0;
        for (const auto& kv : covisibility) {
            maxCovis = std::max(maxCovis, kv.second);
        }

        std::cout << "  Building candidate pairs..." << std::endl;
        for (const auto& kv : covisibility) {
            const auto& pair = kv.first;
            int count = kv.second;

            if (count < config_.minCovisCount)
                continue;

            auto itA = photos.find(pair.first);
            auto itB = photos.find(pair.second);
            if (itA == photos.end() || itB == photos.end())
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

            // Calculate score - adjusted formula
            double overlapScore = (double)cp.covisCount / maxCovis;
            // Use 60° threshold for angle scoring (less aggressive penalty)
            double angleScore = 1.0 - std::min(cp.convergenceAngle / 60.0, 1.0);
            double gsdScore = 1.0 - std::min(std::abs(cp.gsdRatio - 1.0) / 0.5, 1.0);
            double dirScore = std::max(0.0, cp.directionSimilarity);

            cp.score = config_.weight_overlap * overlapScore + config_.weight_angle * angleScore +
                       config_.weight_gsd * gsdScore + config_.weight_direction * dirScore;

            candidates.push_back(cp);
        }

        std::cout << "  Candidate pairs: " << candidates.size() << std::endl;

        // 2. Sort by score
        std::sort(candidates.begin(), candidates.end());

        // 3. Build candidate list for each photo
        std::unordered_map<int, std::vector<const CandidatePair*>> photoCandidates;
        for (const auto& cp : candidates) {
            photoCandidates[cp.photoA].push_back(&cp);
            photoCandidates[cp.photoB].push_back(&cp);
        }

        // 4. Greedy selection + redundancy detection
        std::cout << "  Smart selection..." << std::endl;
        PairGraph selectedGraph;
        std::set<std::pair<int, int>> selectedPairs;
        std::unordered_map<int, int> photoNeighborCount;

        int targetEdges = (int)(photos.size() * config_.targetRedundancyRatio);

        for (const auto& cp : candidates) {
            if ((int)selectedPairs.size() >= targetEdges)
                break;

            int a = cp.photoA;
            int b = cp.photoB;

            // Check neighbor limit
            if (photoNeighborCount[a] >= config_.maxNeighborsPerPhoto &&
                photoNeighborCount[b] >= config_.maxNeighborsPerPhoto) {
                continue;
            }

            // Redundancy detection: skip if already 2-hop connected with enough common neighbors
            if (config_.enableRedundancyCheck && selectedGraph.is2HopConnected(a, b)) {
                int commonNeighbors = selectedGraph.countCommonNeighbors(a, b);
                if (commonNeighbors >= config_.maxCommonNeighborsForRedundancy) {
                    // But still need this edge if one side has too few neighbors
                    if (photoNeighborCount[a] >= config_.minNeighborsPerPhoto &&
                        photoNeighborCount[b] >= config_.minNeighborsPerPhoto) {
                        continue;
                    }
                }
            }

            // Select this edge
            selectedPairs.insert({a, b});
            selectedGraph.addEdge(a, b);
            photoNeighborCount[a]++;
            photoNeighborCount[b]++;
        }

        // 5. Ensure minimum neighbors per photo
        std::cout << "  Ensuring connectivity..." << std::endl;
        for (const auto& kv : photoCandidates) {
            int photoId = kv.first;
            const auto& candidateList = kv.second;

            if (photoNeighborCount[photoId] < config_.minNeighborsPerPhoto) {
                for (const auto* cp : candidateList) {
                    if (photoNeighborCount[photoId] >= config_.minNeighborsPerPhoto)
                        break;

                    int a = std::min(cp->photoA, cp->photoB);
                    int b = std::max(cp->photoA, cp->photoB);

                    if (selectedPairs.count({a, b}) == 0) {
                        selectedPairs.insert({a, b});
                        selectedGraph.addEdge(a, b);
                        photoNeighborCount[cp->photoA]++;
                        photoNeighborCount[cp->photoB]++;
                    }
                }
            }
        }

        // 6. Angle diversity: ensure some pairs with larger angles for geometric diversity
        if (config_.enableAngleDiversity) {
            std::cout << "  Adding angle diversity..." << std::endl;
            std::unordered_map<int, bool> hasLargeAnglePair;

            // Check existing pairs for angle diversity
            for (const auto& cp : candidates) {
                int a = std::min(cp.photoA, cp.photoB);
                int b = std::max(cp.photoA, cp.photoB);
                if (selectedPairs.count({a, b}) > 0 &&
                    cp.convergenceAngle >= config_.angleDiversityThreshold) {
                    hasLargeAnglePair[cp.photoA] = true;
                    hasLargeAnglePair[cp.photoB] = true;
                }
            }

            // Add large angle pairs for photos that don't have any
            for (const auto& kv : photoCandidates) {
                int photoId = kv.first;
                if (hasLargeAnglePair[photoId])
                    continue;

                const auto& candidateList = kv.second;
                for (const auto* cp : candidateList) {
                    if (cp->convergenceAngle >= config_.angleDiversityThreshold &&
                        cp->covisCount >= 50) {  // Require minimum overlap
                        int a = std::min(cp->photoA, cp->photoB);
                        int b = std::max(cp->photoA, cp->photoB);

                        if (selectedPairs.count({a, b}) == 0) {
                            selectedPairs.insert({a, b});
                            selectedGraph.addEdge(a, b);
                            photoNeighborCount[cp->photoA]++;
                            photoNeighborCount[cp->photoB]++;
                            hasLargeAnglePair[cp->photoA] = true;
                            hasLargeAnglePair[cp->photoB] = true;
                            break;
                        }
                    }
                }
            }
        }

        // Convert to result
        std::vector<std::pair<int, int>> result(selectedPairs.begin(), selectedPairs.end());

        // Print statistics
        printStatistics(result, photos, selectedGraph);

        return result;
    }

  private:
    void printStatistics(const std::vector<std::pair<int, int>>& pairs,
                         const std::map<int, Photo>& photos, const PairGraph& graph) {
        std::cout << "\n========== Selection Statistics ==========" << std::endl;
        std::cout << "  Selected pairs: " << pairs.size() << std::endl;
        std::cout << "  Photos involved: " << graph.nodeCount() << std::endl;

        double avgDegree = 2.0 * pairs.size() / graph.nodeCount();
        std::cout << "  Average degree: " << avgDegree << std::endl;

        int mstEdges = (int)graph.nodeCount() - 1;
        double redundancyRatio = (double)pairs.size() / mstEdges;
        std::cout << "  Redundancy ratio: " << redundancyRatio << "x MST" << std::endl;

        double cc = graph.calcClusteringCoefficient();
        std::cout << "  Clustering coefficient: " << cc << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }

    Config config_;
};

// ============================================================================
// Output Writer
// ============================================================================

void writePairs(const std::string& outputPath, const std::vector<std::pair<int, int>>& pairs,
                int totalPhotos) {
    std::ofstream file(outputPath);
    if (!file.is_open()) {
        std::cerr << "Cannot write to file: " << outputPath << std::endl;
        return;
    }

    file << "#Graph of Views" << std::endl;
    file << "#Version 2" << std::endl;
    file << totalPhotos << std::endl;

    for (const auto& p : pairs) {
        file << p.first << " " << p.second << std::endl;
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
    std::cout << "   Pair Selection Algorithm v1.0" << std::endl;
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
    // Parameters adjusted to match reference pair counts
    config.weight_overlap = 0.80;           // Higher overlap priority
    config.weight_angle = 0.05;             // Minimal angle penalty
    config.weight_gsd = 0.10;
    config.weight_direction = 0.05;
    config.minNeighborsPerPhoto = 2;
    config.maxNeighborsPerPhoto = 4;
    config.targetRedundancyRatio = 1.70;    // Higher to match reference counts
    config.enableRedundancyCheck = true;
    config.maxCommonNeighborsForRedundancy = 2;  // Slightly relaxed
    config.enableAngleDiversity = false;
    config.angleDiversityThreshold = 10.0;

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

        // Select pairs
        PairSelector selector(config);
        auto selectedPairs = selector.select(parser.getPhotos(), parser.getCovisibility());

        // Write output
        std::string outputPath = dataDir + "/pairs_generated.txt";
        writePairs(outputPath, selectedPairs, (int)parser.getPhotos().size());

        // Compare with original
        auto originalPairs = loadOriginalPairs(pairsPath);
        if (!originalPairs.empty()) {
            compareResults(selectedPairs, originalPairs);
        }
    }

    std::cout << "\nProcessing complete!" << std::endl;
    return 0;
}
