file(GLOB_RECURSE ALL_SOURCE_FILES src/*.c src/*.h src/*.hpp)

set(ENABLE_CPPCHECK OFF CACHE BOOL "Enable cppcheck")

if(${ENABLE_CPPCHECK})
	add_custom_target(
		    cppcheck 
		    ALL
		    COMMAND cppcheck
		    --enable=warning,performance,portability,information
		    --suppress=missingInclude
		    --suppress=missingIncludeSystem
		    --std=c11
		    --library=std.cfg
		    --template="[{severity}][{id}] {message} {callstack} \(On {file}:{line}\)"
		    --verbose
		    --quiet
		    ${CPPCHECK_FILES}
	)
endif()
