.First.lib <-function (lib, pkg)   {
    library.dynam("MaxChi", pkg, lib)
    library.dynam("Phylpro", pkg, lib)
}
