import re
from typing import Any, List, Dict, Union
from itertools import combinations


def collect_matrix_cols(mat):
    """Convert a matrix into a list of tuples (column-wise)."""
    return [tuple(col) for col in zip(*mat)]


def dict_zip(keys, values):
    """Create a dictionary from keys and values."""
    return dict(zip(keys, values))


def fill_missing_keys(x: dict, keys: list, val=None) -> dict:
    """Add missing keys to a dictionary with a default value."""
    result = x.copy()
    for key in keys:
        if key not in result:
            result[key] = val
    return {k: result.get(k, val) for k in keys}


def filter_list_for_match(items: list, patterns: Union[str, list]) -> list:
    """Return elements of a list matching particular substring(s)."""
    if isinstance(patterns, str):
        patterns = [patterns]

    result = items.copy() if isinstance(items, list) else list(items)
    for pattern in patterns:
        result = [item for item in result if pattern in str(item)]
    return result


def find_first_match_index(pattern: str, items: list) -> int:
    """Return the index of the first match given a pattern."""
    for i, item in enumerate(items):
        if re.search(pattern, str(item)):
            return i
    raise ValueError(f"No match found for pattern: {pattern}")


def interleave(a: list, b: list) -> list:
    """Interleave two lists."""
    result = []
    for i in range(max(len(a), len(b))):
        if i < len(a):
            result.append(a[i])
        if i < len(b):
            result.append(b[i])
    return result


def items_in_a_not_b(a: list, b: list) -> list:
    """Return items in a that are not in b."""
    b_set = set(b)
    return [item for item in a if item not in b_set]


def list2matrix(z: list, ncol: int = 2, byrow: bool = True, fill=None) -> list:
    """Rearrange a list into a matrix (list of lists)."""
    remain = len(z) % ncol
    if remain > 0:
        z = list(z) + [fill] * (ncol - remain)

    if byrow:
        return [z[i:i+ncol] for i in range(0, len(z), ncol)]
    else:
        nrow = len(z) // ncol
        return [[z[j*nrow + i] for j in range(ncol)] for i in range(nrow)]


def matrix2list(mat: list) -> dict:
    """Rearrange a matrix to a named dictionary."""
    result = {}
    for i, row in enumerate(mat):
        for j, val in enumerate(row):
            result[f"{i}-{j}"] = val
    return result


def move_list_items_to_front(items: list, ordered: list) -> list:
    """Rearrange items with ordered items at front."""
    front = [x for x in ordered if x in items]
    back = items_in_a_not_b(items, ordered)
    return front + back


def multiple_replacement(items: Union[str, list], replace_dict: dict) -> Union[str, list]:
    """Replace items using regex patterns from a dictionary."""
    is_string = isinstance(items, str)
    if is_string:
        items = [items]

    result = []
    for item in items:
        text = str(item)
        for pattern, replacement in replace_dict.items():
            text = re.sub(pattern, replacement, text)
        result.append(text)

    return result[0] if is_string else result


def replace_specific_items(items: list, replacer: dict) -> list:
    """Replace specific items in a list using exact matches."""
    return [replacer.get(item, item) for item in items]
