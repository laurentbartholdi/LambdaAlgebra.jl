# functions to print tables

function print_curtis_table(io::IO, Λ::Algebra, totaldegree::Int)
    curtis_string(g) = string(is_lambda(g) ? index(g) : -index(g))
    cache_basis!(Λ,totaldegree)
    for top=1:totaldegree
        for hom=1:totaldegree-top
            first = true
            for μ=hom:-1:0
                basis = get(Λ.basis,(μ=μ,λ=hom-μ,top=top),Nothing[])
                for x=basis
                    if !is_tagger(x)
                        if first
                            print(io, top,",",hom)
                            first = false
                        end
                        print(io," 1(",join((curtis_string(g) for g=x.v)," "),")")
                        if is_taggee(x)
                            print(io, "/",x.tag.second,"(", join((curtis_string(g) for g=x.tag.first)," "),")")
                        end
                    end
                end
            end
            first || println(io)
        end
    end
end
print_curtis_table(Λ::Algebra, totaldegree::Int) = print_curtis_table(stdout,Λ,totaldegree)

function print_grid(io::IO, Λ::Algebra{p}; top_degree = nothing, names = nothing, width = nothing, kwargs...) where p
    if isa(top_degree,Int)
        cache_basis!(Λ,top_degree + top_degree÷(2p-2) + 1)
    else
        top_degree = Λ.total[]
    end
    if names==nothing
        names = top_degree ≤ 30
    end
    if width==nothing
        width = names ? 150 : 60
    end

    backend = get(kwargs,:backend,:text)     
    bot_degree = p==2 ? 1 : 3
    
    data = Matrix{String}(undef,top_degree,top_degree)
    maxhom = 1
    fill!(data,"")
    arrows = Pair{NTuple{3,Int},NTuple{3,Int}}[]
    for top=bot_degree:top_degree
        for hom=1:top
            if hom+top > Λ.total[]
                data[top,hom] = "?"
                continue
            end
            entries = String[]
            for μ=0:hom
                basis = get(Λ.basis,(μ=μ,λ=hom-μ,top=top),Nothing[])
                i = 0
                for x=basis
                    i += 1
                    xstring = string(names ? x.v : μ)
                    
                    if is_tagger(x) && backend==:tikz
                        rangedegree = degree(Λ,x.tag.first)
                        j = 0
                        for y=Λ.basis[rangedegree] # find position in range
                            j += 1
                            y.v==x.tag.first && break
                        end
                        push!(arrows,(top,hom,i)=>(top-1,rangedegree.λ+rangedegree.μ,j))
                    end
                    if is_tagger(x) && backend==:html
                        rangedegree = degree(Λ,x.tag.first)
                        shift = rangedegree.λ+rangedegree.μ-(basis.degree.λ+basis.degree.μ)
                        x1,x2 = 5, 5+width*shift
                        y1,y2 = 12,28
                        arrow = """<svg width="$x2" style="position: absolute;
    width: 0;
    height: 0;
    top: 0;
    left: 0;">
  <defs>
    <marker id="redhead" orient="auto" markerWidth="6" markerHeight="4" refX="6" refY="2" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,4 L6,2 Z" fill="red"/>
    </marker>
  </defs>
  <line x1="$x1" y1="$y1" x2="$x2" y2="$y2" style="stroke:rgb(255,0,0);stroke-width:1.5" marker-end="url(#redhead)"/>
</svg>"""
                        push!(entries, arrow*xstring)
                    elseif is_alive(x) && backend==:html
                        push!(entries, """<span style="background-color:chartreuse;">$xstring</span>""")
                    elseif is_alive(x) && backend==:tikz
                        push!(entries, """\\boxed{$xstring}""")
                    else
                        push!(entries, xstring)
                    end
                end
            end
            if !isempty(entries)
                maxhom = max(hom,maxhom)
                data[top,hom] = join(entries,",")
            end
        end
    end
    
    if backend==:latex
        kwargs = (highlighters = [LatexHighlighter((_, i, _)->iszero((top_degree+1-i)%5), ["textbf"])],
                  kwargs...)
    elseif backend==:html
        kwargs = (maximum_column_width = string(width,"px"),
                  allow_html_in_cells = true,
                  style = HtmlTableStyle(first_line_column_label = ["width" => string(width,"px")]),
                  kwargs...)
    elseif backend==:markdown
        kwargs = (highlighters = [MarkdownHighlighter((_, i, _)->iszero((top_degree+1-i)%5), MarkdownStyle(bold=true))],
                  kwargs...)
    else 
        kwargs = (display_size = (-1,-1),
                  
                  highlighters = [TextHighlighter((_, i, _)->iszero((top_degree+1-i)%5), crayon"fg:black bold bg:light_gray")],
                  #table_format = TextTableFormat(horizontal_lines_at_data_rows=collect(mod1(top_degree,5):5:top_degree)),
                  kwargs...)
    end

    if backend==:tikz
        println(io,"""\\documentclass{standalone}
\\usepackage{tikz}
\\usetikzlibrary{matrix,positioning,calc}
\\usepackage{unicode-math,amsmath}
\\begin{document}
\\tikzset{toprule/.style={%
        execute at end cell={%
            \\draw [line cap=rect,#1] (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.north west) -- (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.north east);%
        }
    },
    bottomrule/.style={%
        execute at end cell={%
            \\draw [line cap=rect,#1] (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.south west) -- (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.south east);%
        }
    }
}
\\begin{tikzpicture}[node distance=0.5em,>=stealth]
\\matrix[row sep=2ex,column sep=1em] {""")
        for i=1:maxhom
            print(io," & \\node{\\textbf{",i,"}};")
        end
        println(io, "\\\\")
        for i=top_degree:-1:bot_degree
            print(io,"\\node[minimum height=3ex,rectangle] (c",i,"-begin) {\\textbf{",i,"}};")
            for j=1:maxhom
                s = split(data[i,j],",")
                print(io,"&")
                placement = ""
                count = 1
                for a=s
                    print(io," \\node[inner sep=0pt,outer sep=0pt,",placement,"] (c",i,"-",j,"-",count,") {\\(",a,"\\)};")
                    placement = "right=of c$i-$j-$count"
                    count += 1
                end
                println(io,"\\node[",placement,"] (c",i,"-",j,"-end) {};")
            end
            println(io,"\\\\")
        end
        println(io,"};")
        for (a,b)=arrows
            println(io,"\\draw[red,->] (c",a[1],"-",a[2],"-",a[3],".south) -- (c",b[1],"-",b[2],"-",b[3],".north);")
        end
        for i=mod1(top_degree,5):5:top_degree
            println(io,"\\draw (c$i-begin.north west) -- (c$i-$maxhom-end.north east |- c$i-begin.north west);")
        end
        println(io,"""\\end{tikzpicture}
\\end{document}""")
        return
    end
    
    pretty_table(io, data[top_degree:-1:bot_degree,1:maxhom];
                 column_labels=1:maxhom, row_labels=top_degree:-1:bot_degree,
                 title = "$Λ",
                 kwargs...)
end
print_grid(Λ::Algebra; kwargs...) = print_grid(stdout, Λ; kwargs...)

function pdf_grid(name::AbstractString, Λ::Algebra{p}; top_degree = nothing, names = nothing) where p
    mktempdir() do dir
        cd(dir) do
            tex = tempname(dir,suffix=".tex")
            open(tex,"w") do f
                print_grid(f, Λ, top_degree=top_degree, names=names, backend=:tikz)
            end
            out = IOBuffer()
            if !success(pipeline(`lualatex $tex`,stdout=out)) # lualatex because of memory usage
                @error String(take!(out))
            end
            mv(tex[1:end-3]*"pdf",name,force=true)
        end
    end
end

function print_stem(io::IO,Λ::Vector{Algebra{p,Names}},k::Int,full=true;odd=false) where {p,Names}
    vars = Monomial{p,Names}[]
    step = 1+Int(odd)
    if full
        for i=1:step:length(Λ)
            print(io,rpad(i,5))
        end
        println(io)
    end
    for L=Λ[1:step:length(Λ)]
        s = ""
        for b=values(L.basis)
            b.degree.top==k || continue
            for x=b
                if is_alive(x)
                    if x.v∉vars
                        push!(vars,x.v)
                    end
                    i = findfirst(==(x.v),vars)
                    if full
                        s *= 'a'+i-1
                    else
                        s *= '*'
                    end
                end
            end
        end
        print(io,rpad(s,5))
    end
    println(io)
    if full
        for i=1:length(vars)
            println(io,'a'+i-1,": ",vars[i]," ",degree(Λ[1],vars[i]))
        end
    end
end
print_stem(Λ::Vector{Algebra{p,Names}},k::Int;kwargs...) where {p,Names} = print_stem(stdout,Λ,k;kwargs...)

function print_stem(io::IO,Λ::Vector{Algebra{p,Names}},k::AbstractRange;kwargs...) where {p,Names}
    print(io,rpad("*",5))
    step = Int(get(kwargs,:odd,false))+1
    for i=1:step:length(Λ)
        print(io,rpad(i,5))
    end
    println(io)
    for k=k
        print(io,rpad(k,5))
        print_stem(io,Λ,k,false;kwargs...)
    end
end
    
print_stem(Λ::Vector{Algebra{p,Names}},k::AbstractRange;kwargs...) where {p,Names} = print_stem(stdout,Λ,k;kwargs...)
