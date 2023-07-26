import Literate

corefiles = sort(filter(f -> startswith(f, r"\d\d"), readdir()))

vignetteconfig = Dict(
    "repo_root_url" => "https://github.com/PoisotLab/Omnomnomivores",
    "codefence" => Pair("````julia", "````"),
    "flavor" => Literate.FranklinFlavor(),
    "credit" => false
)

nbconfig = Dict(
    "execute" => false
)

for corefile in corefiles
    Literate.markdown(corefile, "vignettes"; config=vignetteconfig)
    Literate.notebook(corefile, "vignettes"; config=nbconfig)
end