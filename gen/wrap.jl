using Clang
using Clang.Generators

headers = ["../include/soil.h", "../include/beps.h"]
args = get_default_args()

options = load_options(joinpath(@__DIR__, "generator.toml"))

ctx = create_context(headers, args)
build!(ctx)
