# Package
version       = "0.01"
author        = "Brent Pedersen"
description   = "fast pileup"
license       = "MIT"


# Dependencies
skipDirs = @["tests", "src/hilepkg"]

requires "hts >= 0.2.13"
srcDir = "src"
installExt = @["nim"]


