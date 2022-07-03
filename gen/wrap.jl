using Clang
using Clang.Generators

headers = ["../include/soil.h", "../include/beps.h"]
args = get_default_args()



ctx = create_context(headers, args)
build!(ctx)
