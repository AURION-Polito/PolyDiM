file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.h *.hpp)

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "Enable clangformat")

set(CLANGFORMAT_STYLE "{BasedOnStyle: Microsoft, Language: Cpp, Standard: Cpp11, PointerAlignment: Left, DerivePointerAlignment: false, TabWidth: 4, IndentWidth: 4, UseTab: Always, AlignAfterOpenBracket: true, BreakBeforeBraces: Custom, SortUsingDeclarations: true, SortIncludes: true}")

if(${ENABLE_CLANGFORMAT})
	add_custom_target(
		clangformat
		ALL
		COMMAND clang-format
		-style=${CLANGFORMAT_STYLE}
		-i
		${CLANGFORMAT_FILES}
	)
endif()
