# NEXUS trees file parser (e.g. as returned by RevBayes `writeNexus`)
"""
    readnex(fname)

Read a Nexus file. This is hacky and not well-tested, only for
RevBayes output files... Let's say it is 'unoffical' as of yet.
"""
function readnex(fname)
    content = readuntil(fname, "End;")
    treelines = split(split(content, "Begin")[2], "\n")[2:end]
    treelines = filter(x->x!="", treelines)
    map(treelines) do treeline
        tree = strip(join(split(treeline, "=")[2:end]))
        readnw(replace(tree, r"\[.+?\]"=>""))
    end
end
