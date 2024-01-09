# List of Examples

```@eval
# if we ever need something fancier than @contents, e.g., including
# tags, the following provides a starting point:
using Markdown

lines = []
for (title, index_md) in Main.EXAMPLES
    folder = splitpath(index_md)[2]
    push!(lines, "* [$title]($folder/)")
end
Markdown.parse(join(lines, "\n"))
nothing # XXX
```

```@contents
Pages = [index_md for (title, index_md) in EXAMPLES]
Depth = 1
```
