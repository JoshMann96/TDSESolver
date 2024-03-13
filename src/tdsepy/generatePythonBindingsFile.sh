#!/usr/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_cpp_file> <output_python_file>"
    exit 1
fi

CPP_FUNC_REGEX="^\s*(\w+\*?)\s+(\w+)\s*\(\s*([^)]*)\s*\)"  # functionName(argType0 argName0 argType1 argName1...)
# ^(\w+\*?) - outputType
# \s+    - 1+ whitespace
# (\w+)  - functionName
# \s*\(\s*  - 0+ white spaces around parenthesis
# ([^)]*) - parameters
# \s*\)\s*  - 0+ white spaces around right parenthesis

ARG_REGEX="(\w+\*?)\s+(\w+)"


function PRINT_TYPESET_START {
    echo -n "lib.$1.argtypes = [" >&3
}

function PRINT_TYPESET_END {
    echo "]" >&3
}

function PRINT_DEF_START {
    echo -n "def $1(" >&3
}

function PRINT_DEF_END {
    echo "):" >&3
}

function PRINT_CALL_START {
    if [[ "$2" == "void*" ]]; then
        echo -n "    return ctypes.c_void_p(lib.$1(" >&3
    else
        echo -n "    return lib.$1(" >&3
    fi
}

function PRINT_CALL_END {
    if [[ "$1" == "void*" ]]; then
        echo "))" >&3
    else
        echo ")" >&3
    fi
}

function PRINT_CTYPE {
    case "$1" in
        "int")
            echo -n "ctypes.c_int" >&3
            ;;
        "double")
            echo -n "ctypes.c_double" >&3
            ;;
        "void")
            echo -n "ctypes.c_void" >&3
            ;;
        "void*")
            echo -n "ctypes.c_void_p" >&3
            ;;
        "string")
            echo -n "ctypes.c_wchar_p" >&3
            ;;
    esac
}

function PRINT_PTYPE {
    case "$1" in
        "int")
            echo -n "int" >&3
            ;;
        "double")
            echo -n "float" >&3
            ;;
        "void*")
            echo -n "ctypes.c_void_p" >&3
            ;;
        "string")
            echo -n "str" >&3
            ;;
    esac
}

function HANDLE_TYPE_SETTING {
    MSG=$1
    I=$2

    if [[ $MSG =~ $ARG_REGEX ]] ; then
        DTYPE="${BASH_REMATCH[1]}"
        ARGST="${BASH_REMATCH[2]}"

        if [[ $I != 0 ]] ; then
            echo -n ", " >&3
        fi
        PRINT_CTYPE "$DTYPE"

        # Remove the first regex match and try again
        HANDLE_TYPE_SETTING "${MSG/${BASH_REMATCH[0]}/}" ${I+1}
    fi
}

function HANDLE_ARG_MATCHES {
    MSG=$1
    I=$2

    if [[ $MSG =~ $ARG_REGEX ]] ; then
        DTYPE="${BASH_REMATCH[1]}"
        ARGST="${BASH_REMATCH[2]}"


        if [[ $I != 0 ]] ; then
            echo -n ", " >&3
        fi
        echo -n "$ARGST:" >&3
        PRINT_PTYPE "$DTYPE"


        # Remove the first regex match and try again
        HANDLE_ARG_MATCHES "${MSG/${BASH_REMATCH[0]}/}" ${I+1}
    fi
}

function HANDLE_CTYPE_MATCHES {
    MSG=$1

    if [[ $MSG =~ $ARG_REGEX ]] ; then
        DTYPE="${BASH_REMATCH[1]}"
        ARGST="${BASH_REMATCH[2]}"

        if [[ "$DTYPE" != "void*" ]]; then
            echo -n "    ${ARGST} = " >&3
            PRINT_CTYPE "$DTYPE"
            echo "(${ARGST})" >&3
        fi

        # Remove the first regex match and try again
        HANDLE_CTYPE_MATCHES "${MSG/${BASH_REMATCH[0]}/}"
    fi
}

function HANDLE_CALL_ARGS {
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
        echo -n "lib.${FUNCN}.restype = " >&3
        PRINT_CTYPE "$DTYPE"
        echo "" >&3
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

# Open the output Python file for writing
exec 3> "$2"

cat _PyBindingsPreamble >&3

IN_PYDOC=0

PYDOC_BUF=$(mktemp)

#Go through each line in file
while IFS= read -r line; do
    # Starting/ending docstring?
    if [[ $line =~ \s*\"{3} ]]; then
        IN_PYDOC=$((!IN_PYDOC))
        echo "${line}" >> "$PYDOC_BUF"
    # Currently doing docstring?
    elif [[ ${IN_PYDOC} == 1 ]]; then
        echo "${line}" >> "$PYDOC_BUF"
    # FUNCTION DEFINITION
    elif [[ $line =~ $CPP_FUNC_REGEX ]]; then
        # Extract components from the matched line
        outputType="${BASH_REMATCH[1]}"       # outputType
        functionName="${BASH_REMATCH[2]}"     # functionName
        functionArgs="${BASH_REMATCH[3]}"             # argType0 argName0 argType1 argName1...

        ### RUN ONLY IF IS FILE DEFINITION ###
        if [[ "${outputType}" != "return" ]]; then
            echo "" >&3

            ### CTYPES SETTING ###
            PRINT_TYPESET "$functionName" "$functionArgs"
            PRINT_RESTYPESET "$functionName" "$outputType"

            echo "" >&3

            ### FUNCTION DEFINITION ###
            PRINT_DEF_START "$functionName"
            PRINT_DEF_ARGS "$functionArgs"
            PRINT_DEF_END

            ### PRINT DOCSTRING IN BUFFER ###
            cat "$PYDOC_BUF" >&3
            #clear docstring buffer
            > "$PYDOC_BUF"

            ### CTYPE CONVERSIONS ###
            PRINT_CTYPE_CONV "$functionArgs"

            ### CALL FUNCTION ###
            PRINT_CALL_START "$functionName" "$outputType"
            PRINT_CALL_ARGS "$functionArgs"
            PRINT_CALL_END "$outputType"
        fi
    fi
done < "$1"

# Close the output Python file
exec 3>&-

rm "$PYDOC_BUF"

exit