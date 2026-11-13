-- === Wells + controls added by build_deck.py ===
WELSPECS
-- Name    Group   I   J   RefDepth    Phase
   'I01'    'INJ'   6     4     1*      WATER /
   'I03'    'INJ'   17     13     1*      WATER /
   'I04'    'INJ'   7     20     1*      WATER /
   'I02'    'INJ'   23     28     1*      WATER /
   'P06'    'PROD'   3     27     1*      OIL /
   'P01'    'PROD'   26     3     1*      OIL /
   'P02'    'PROD'   23     11     1*      OIL /
   'P03'    'PROD'   17     18     1*      OIL /
   'P05'    'PROD'   19     2     1*      OIL /
   'P04'    'PROD'   11     29     1*      OIL /
/

COMPDAT
-- Name    I   J   K1  K2  State   SatTab  TranFactor  WBDIA
  'I01'    6     4     1    2    OPEN    1*  1* 0.3048 /
  'I01'    6     4     4    7    OPEN    1*  1* 0.3048 /
  'I03'    17     13     1    6    OPEN    1*  1* 0.3048 /
  'I04'    7     20     1    2    OPEN    1*  1* 0.3048 /
  'I04'    7     20     4    7    OPEN    1*  1* 0.3048 /
  'I02'    23     28     4    7    OPEN    1*  1* 0.3048 /
  'P06'    3     27     4    7    OPEN    1*  1* 0.3048 /
  'P06'    3     27     1    1    OPEN    1*  1* 0.3048 /
  'P01'    26     3     4    7    OPEN    1*  1* 0.3048 /
  'P02'    23     11     1    4    OPEN    1*  1* 0.3048 /
  'P03'    17     18     1    3    OPEN    1*  1* 0.3048 /
  'P03'    17     18     5    7    OPEN    1*  1* 0.3048 /
  'P05'    19     2     1    2    OPEN    1*  1* 0.3048 /
  'P05'    19     2     4    6    OPEN    1*  1* 0.3048 /
  'P04'    11     29     1    3    OPEN    1*  1* 0.3048 /
  'P04'    11     29     4    7    OPEN    1*  1* 0.3048 /
/

WCONPROD
-- Name   OpenShut  Ctrl  ORAT    WRAT    GRAT    LRAT    RESV    BHP
   'P01'   OPEN    ORAT    200 4*  50 /
   'P02'   OPEN    ORAT    200 4*  50 /
   'P03'   OPEN    ORAT    200 4*  50 /
   'P04'   OPEN    ORAT    200 4*  50 /
   'P05'   OPEN    ORAT    200 4*  50 /
   'P06'   OPEN    ORAT    200 4*  50 /
/

WCONINJE
-- Name    Type    OpenShut    Ctrl    RATE    RESV    BHP
   'I01'   WATER   OPEN    RATE    300 1*  400 /
   'I02'   WATER   OPEN    RATE    300 1*  400 /
   'I03'   WATER   OPEN    RATE    300 1*  400 /
   'I04'   WATER   OPEN    RATE    300 1*  400 /
/

TSTEP
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31
/

-- === original Petrel schedule content below ===
-- === Wells + controls (CLEAN rebuild, single copy) ===
WELSPECS
-- Name    Group   I   J   RefDepth    Phase
   'I01'    'INJ'   6     4     1*      WATER /
   'I03'    'INJ'   17     13     1*      WATER /
   'I04'    'INJ'   7     20     1*      WATER /
   'I02'    'INJ'   23     28     1*      WATER /
   'P06'    'PROD'   3     27     1*      OIL /
   'P01'    'PROD'   26     3     1*      OIL /
   'P02'    'PROD'   23     11     1*      OIL /
   'P03'    'PROD'   17     18     1*      OIL /
   'P05'    'PROD'   19     2     1*      OIL /
   'P04'    'PROD'   11     29     1*      OIL /
/

COMPDAT
-- Name    I   J   K1  K2  State   SatTab  TranFactor  WBDIA
  'I01'    6     4     1    2    OPEN    1*  1* 0.3048 /
  'I01'    6     4     4    7    OPEN    1*  1* 0.3048 /
  'I03'    17     13     1    6    OPEN    1*  1* 0.3048 /
  'I04'    7     20     1    2    OPEN    1*  1* 0.3048 /
  'I04'    7     20     4    7    OPEN    1*  1* 0.3048 /
  'I02'    23     28     4    7    OPEN    1*  1* 0.3048 /
  'P06'    3     27     4    7    OPEN    1*  1* 0.3048 /
  'P06'    3     27     1    1    OPEN    1*  1* 0.3048 /
  'P01'    26     3     4    7    OPEN    1*  1* 0.3048 /
  'P02'    23     11     1    4    OPEN    1*  1* 0.3048 /
  'P03'    17     18     1    3    OPEN    1*  1* 0.3048 /
  'P03'    17     18     5    7    OPEN    1*  1* 0.3048 /
  'P05'    19     2     1    2    OPEN    1*  1* 0.3048 /
  'P05'    19     2     4    6    OPEN    1*  1* 0.3048 /
  'P04'    11     29     1    3    OPEN    1*  1* 0.3048 /
  'P04'    11     29     4    7    OPEN    1*  1* 0.3048 /
/

WCONPROD
-- Name   OpenShut  Ctrl  ORAT    WRAT    GRAT    LRAT    RESV    BHP
   'P01'   OPEN    ORAT    200 4*  50 /
   'P02'   OPEN    ORAT    200 4*  50 /
   'P03'   OPEN    ORAT    200 4*  50 /
   'P04'   OPEN    ORAT    200 4*  50 /
   'P05'   OPEN    ORAT    200 4*  50 /
   'P06'   OPEN    ORAT    200 4*  50 /
/

WCONINJE
-- Name    Type    OpenShut    Ctrl    RATE    RESV    BHP
   'I01'   WATER   OPEN    RATE    300 1*  400 /
   'I02'   WATER   OPEN    RATE    300 1*  400 /
   'I03'   WATER   OPEN    RATE    300 1*  400 /
   'I04'   WATER   OPEN    RATE    300 1*  400 /
/

TSTEP
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31 91.31 91.31 91.31 91.31
 91.31 91.31 91.31 91.31
/

-- === original Petrel schedule content below ===
-- Generated [
-- Format      : ECLIPSE keywords (ASCII)
-- Exported by : Petrel 2018.2 Schlumberger
-- User name   : cathe
-- Date        : Monday, July 13 2026 00:02:44
-- Project     : try1.pet
-- Generated ]

RPTSCHED                               -- Generated : Petrel
  FIP WELLS /

RPTRST                                 -- Generated : Petrel
  BASIC=3 FLOWS FREQ /

