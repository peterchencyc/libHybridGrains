#ifndef COMPILE_DEFINITIONS_H
#define COMPILE_DEFINITIONS_H

namespace CompileDefinitions {

// SHA1 hash of the Git commit for this build
constexpr char GitSHA1[] = "GITDIR-NOTFOUND";

// Build mode (Debug, Release, etc)
constexpr char BuildMode[] = "Coverage";

// C compiler used for this build
constexpr char CCompiler[] = "/usr/local/bin/gcc-9";

// C++ compiler used for this build
constexpr char CXXCompiler[] = "/usr/local/bin/g++-9";

// Date and time
constexpr char BuildDateTime[] = "2020-01-03 17:11:28";

} // namespace CompileDefinitions

#endif
