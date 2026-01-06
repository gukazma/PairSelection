set(TARGET_NAME QTDemoApp)
find_package(Qt6 REQUIRED Core Gui Widgets)  
qt_standard_project_setup()


file(GLOB_RECURSE CPP_HEAD_FILES  ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB_RECURSE CPP_SOUCE_FILES  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE UI_SOUCE_FILES  ${CMAKE_CURRENT_SOURCE_DIR}/*.ui)
set(QRC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/QResources/resources.qrc)
set(APPRC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/QResources/appicon.rc)
set(Translatefile ${CMAKE_CURRENT_SOURCE_DIR}/QResources/Translations/zh_CN.ts)

source_group(TREE ${CMAKE_CURRENT_LIST_DIR} FILES ${CPP_HEAD_FILES} ${CPP_SOUCE_FILES})
# file(GLOB_RECURSE SHADER_FILES_ABS  ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*)
# set(SHADER_FILES_REL "")
# foreach(SHADER_FILE_ABS ${SHADER_FILES_ABS})
#     file(RELATIVE_PATH SHADER_FILE_REL ${CMAKE_CURRENT_SOURCE_DIR} ${SHADER_FILE_ABS})
#     list(APPEND SHADER_FILES_REL ${SHADER_FILE_REL})
# endforeach()

qt_add_executable(${TARGET_NAME}  ${CPP_HEAD_FILES} ${CPP_SOUCE_FILES} ${UI_SOUCE_FILES}  ${QRC_FILE} ${APPRC_FILE} ${Translatefile})
set_target_properties(${TARGET_NAME} PROPERTIES FOLDER Apps)


target_link_libraries(${TARGET_NAME}
PUBLIC
Qt6::Widgets
Qt6::Gui
)

# qt_add_resources(${TARGET_NAME} "images"
#     PREFIX "/images"
#     FILES ./QResources/mgukazma.jpg
# )

# qt_add_shaders(${TARGET_NAME} "shaders"
# GLSL "310es, 430"
# PREFIX
# "/"
# FILES
# ${SHADER_FILES_REL}
# )

qt_generate_deploy_app_script(
    TARGET ${TARGET_NAME}
    OUTPUT_SCRIPT deploy_script
    NO_UNSUPPORTED_PLATFORM_ERROR
)
install(SCRIPT ${deploy_script})