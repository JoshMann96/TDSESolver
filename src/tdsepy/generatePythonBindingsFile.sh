#!/usr/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_cpp_file> <output_python_file>"
    exit 1
fi

CPP_REGEX="^\s*(\w+)\s+(\w+)\s*\(\s*([^)]*)\s*\)"  # functionName(argType0 argName0 argType1 argName1...)
# ^(\w+) - outputType
# \s+    - 1+ whitespace
# (\w+)  - functionName
# \s*\(\s*  - 0+ white spaces around parenthesis
# ([^)]*) - parameters
# \s*\)\s*  - 0+ white spaces around right parenthesis

PY_DEF_TEMP="def %s(%s):"

# Open the output Python file for writing
exec 3> "$2"

function PRINT_TYPESET_START {
    echo -n "lib.$1.argtypes = [" >&3
}

function PRINT_TYPESET_END {
    echo "]" >&3
}

function PRINT_DEF_START {
    echo -n "def $1(" >&3
}

function PRINT_ARGS {
    echo -n "$1:$2" >&3
}

function PRINT_DEF_END {
    echo "):" >&3
}

function PRINT_CALL_START {
    echo -n "    return lib.$1(" >&3
}

function PRINT_CALL_END {
    echo ")" >&3
}

function HANDLE_TYPE_SETTING {
    ARG_REGEX="(\w+)\s+(\w+)"
    MSG=$1
    I=$2

    if [[ $MSG =~ $ARG_REGEX ]] ; then
        DTYPE="${BASH_REMATCH[1]}"
        ARGST="${BASH_REMATCH[2]}"


        if [[ $I != 0 ]] ; then
            echo -n ", " >&3
        fi
        echo -n "ctypes.c_${DTYPE}" >&3


        # Remove the first regex match and try again
        HANDLE_TYPE_SETTING "${MSG/${BASH_REMATCH[0]}/}" ${I+1}
    fi
}

function HANDLE_ARG_MATCHES {
    ARG_REGEX="(\w+)\s+(\w+)"
    MSG=$1
    I=$2

    if [[ $MSG =~ $ARG_REGEX ]] ; then
        DTYPE="${BASH_REMATCH[1]}"
        ARGST="${BASH_REMATCH[2]}"


        if [[ $I != 0 ]] ; then
            echo -n ", " >&3
        fi
        echo -n "$ARGST" >&3


        # Remove the first regex match and try again
        HANDLE_ARG_MATCHES "${MSG/${BASH_REMATCH[0]}/}" ${I+1}
    fi
}

function HANDLE_CTYPE_MATCHES {
    ARG_REGEX="(\w+)\s+(\w+)"
    MSG=$1

    if [[ $MSG =~ $ARG_REGEX ]] ; then
        DTYPE="${BASH_REMATCH[1]}"
        ARGST="${BASH_REMATCH[2]}"


        echo "    $ARGST = ctypes.c_$DTYPE($ARGST)" >&3


        # Remove the first regex match and try again
        HANDLE_CTYPE_MATCHES "${MSG/${BASH_REMATCH[0]}/}"
    fi
}

function HANDLE_CALL_ARGS {
    ARG_REGEX="(\w+)\s+(\w+)"
    MSG=$1
    I=$2

    if [[ $MSG =~ $ARG_REGEX ]] ; then
        DTYPE="${BASH_REMATCH[1]}"
        ARGST="${BASH_REMATCH[2]}"


        if [[ $I != 0 ]] ; then
            echo -n ", " >&3
        fi
        echo -n "$ARGST" >&3


        # Remove the first regex match and try again
        HANDLE_CALL_ARGS "${MSG/${BASH_REMATCH[0]}/}" ${I+1}
    fi
}

function PRINT_TYPESET {
    PRINT_TYPESET_START "$1"
    HANDLE_TYPE_SETTING "$2" 0
    PRINT_TYPESET_END
}

function PRINT_RESTYPESET {
    FUNCN=$1
    DTYPE=$2
    if [[ "$DTYPE" == "void" ]]; then
        echo "lib.${FUNCN}.restype = None" >&3
    else
        echo "lib.${FUNCN}.restype = ctypes.c_${DTYPE}" >&3
    fi
}

function PRINT_DEF_ARGS {
    HANDLE_ARG_MATCHES "$1" 0
}

function PRINT_CTYPE_CONV {
    HANDLE_CTYPE_MATCHES "$1"
}

function PRINT_CALL_ARGS {
    HANDLE_CALL_ARGS "$1" 0
}

############ BEGIN PYTHON ENCODING #############
#imports
echo "import ctypes" >&3
echo "# Load shared library" >&3
echo "lib = ctypes.CDLL('./libTDSEpy.so')" >&3
echo "">&3


#definitions
while IFS= read -r line; do
    # Check if the line matches the function definition pattern
    if [[ $line =~ $CPP_REGEX ]]; then
        # Extract components from the matched line
        outputType="${BASH_REMATCH[1]}"       # outputType
        functionName="${BASH_REMATCH[2]}"     # functionName
        functionArgs="${BASH_REMATCH[3]}"             # argType0 argName0 argType1 argName1...

        

        echo "" >&3

        ### CTYPES SETTING ###
        PRINT_TYPESET "$functionName" "$functionArgs"
        PRINT_RESTYPESET "$functionName" "$outputType"

        echo "" >&3

        ### FUNCTION DEFINITION ###
        PRINT_DEF_START "$functionName"
        PRINT_DEF_ARGS "$functionArgs"
        PRINT_DEF_END

        PRINT_CTYPE_CONV "$functionArgs"

        ### CALL FUNCTION ###
        PRINT_CALL_START "$functionName"
        PRINT_CALL_ARGS "$functionArgs"
        PRINT_CALL_END
    fi
done < "$1"

# Close the output Python file
exec 3>&-

exit