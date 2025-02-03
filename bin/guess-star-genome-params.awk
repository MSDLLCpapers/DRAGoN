#!/usr/bin/awk -f

# Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
# This file is part of DRAGoN.
#
# This source code is licensed under the MIT License found in the
# LICENSE file in the root directory of this source tree.

/^>/ {nchrom++}
!/^>/ {glen += length($0)}
END {
    print "#glen = " glen ", nchrom = " nchrom
    genomeSAindexNbases = int(log(glen) / log(2) - 1)
    genomeChrBinNbits = int(log(glen / nchrom > 150 ? glen / nchrom : 150) / log(2))
    print "export genomeSAindexNbases=" (genomeSAindexNbases < 14 ? genomeSAindexNbases : 14)
    print "export genomeChrBinNbits=" (genomeChrBinNbits < 18 ? genomeChrBinNbits : 18)
}
