module ClusteringExt

using NewickTree, Clustering
using NewickTree: setdistance!

function Base.convert(::Type{Node}, hc::Hclust)
   nodes = [Node(i, n=string(n), d=0.) for (i,n) in zip(hc.order,hc.labels)]
   n = length(nodes)
   idfun(x) = x > 0 ? x + n : abs(x)
   for i=1:size(hc.merges, 1)
       nid = n + i
       j, k = idfun.(hc.merges[i,:])
       a = nodes[j]
       b = nodes[k]
       h = hc.heights[i]
       newnode = Node(nid, n="$nid", d=h)
       setdistance!(a, h-distance(a))
       setdistance!(b, h-distance(b))
       push!(newnode, a)
       push!(newnode, b)
       push!(nodes, newnode)
   end
   setdistance!(nodes[end], 0.)
   return nodes[end]
end

end
