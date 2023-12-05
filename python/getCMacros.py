from typing import Dict
import re

def getCMacros(fileName: str) -> Dict:
    with open(fileName,"r") as f:
        lines = f.readlines()

    macros = {}
    for line in lines:
        if line.startswith("#define"):
            regex = r'\(.+?\)|\S+'
            result = re.findall(regex, line)
            if len(result) > 2:
                macros[result[1]] = result[2]
    return macros