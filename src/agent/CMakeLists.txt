# ---- stimlab_agent : LLM agent loop (providers, tools, tool-calling loop) ----
# Builds under BOTH presets. The Anthropic provider's networking is guarded by
# STIMLAB_HAVE_SCIENCE; without the science feature it compiles to a stub and the
# agent falls back to the offline MockProvider. So the default `windows` preset
# stays curl-free and fast.
add_library(stimlab_agent STATIC
  Tools.cpp
  SystemPrompt.cpp
  MockProvider.cpp
  AnthropicProvider.cpp
  WebTools.cpp
  Agent.cpp
)

target_include_directories(stimlab_agent PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_agent PUBLIC
  stimlab_contracts
  stimlab_data
  stimlab_core
  nlohmann_json::nlohmann_json
)

# libcurl only in the science build; all curl code is #ifdef STIMLAB_HAVE_SCIENCE.
if(STIMLAB_ENABLE_SCIENCE)
  find_package(CURL REQUIRED)
  target_link_libraries(stimlab_agent PUBLIC CURL::libcurl)
endif()

target_compile_features(stimlab_agent PUBLIC cxx_std_20)
