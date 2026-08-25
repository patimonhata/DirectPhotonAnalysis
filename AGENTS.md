## File editing

- `apply_patch` is unavailable in this environment.
- Do not invoke or probe for `apply_patch`.
- Use `git apply` directly for patch-based edits.
- If `git apply` is unsuitable, use another available editing method without retrying `apply_patch`.

## C++ formatting

- Prefer vertically compact C++ code. Do not wrap lines merely to avoid horizontal scrolling.
- Treat 160 columns as a soft limit, not an absolute limit. A clear, indivisible expression may exceed it moderately.
- Keep a declaration and its initializer on one line when the complete line is reasonably readable.
- Do not split a single file path, URL, object name, branch name, or similar atomic string into adjacent string literals solely because of line length.
- Keep short function calls, argument lists, and boolean expressions on one line when they remain easy to understand.
- Introduce line breaks when they expose meaningful structure, including:
  - long stream output expressions;
  - complex ternary or boolean expressions;
  - long struct/member initializer lists;
  - calls or declarations with many logically distinct arguments;
  - lists where one item per line improves comparison or review.
- When wrapping is necessary, wrap at semantic boundaries and use as few lines as practical.
- Preserve the surrounding style and do not reformat unrelated code.