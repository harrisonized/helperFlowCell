import re


def title_to_snake_case(text: str) -> str:
    """Convert 'Column Title' to 'column_title'."""
    return '_'.join(text.lower().split())


def substr_right(x: str, n: int) -> str:
    """Extract last n characters from a string."""
    return x[-n:] if len(x) >= n else x


def txt_strip(x: str, chars: str = ' ') -> str:
    """Remove leading and trailing characters."""
    char_list = list(set(chars))

    for char in char_list:
        # escape special regex characters
        if char in '()[]{}.*+?^$|\\':
            char = '\\' + char

        x = re.sub(f'^{char}+', '', x)
        x = re.sub(f'{char}+$', '', x)

    return x
