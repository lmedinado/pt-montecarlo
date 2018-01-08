#pragma once

const unsigned int DEFAULT = 0;

const unsigned int NO_START   = 0x0;
const unsigned int COLD_START = 0x1 << 0;
const unsigned int HOT_START  = 0x1 << 1;
const unsigned int START_MASK = COLD_START | HOT_START;

const unsigned int RANDOM_WALK    = 0x0;
const unsigned int UNIFORM_SEARCH = 0x1 << 2;
const unsigned int SEARCH_MASK    = RANDOM_WALK | UNIFORM_SEARCH;

const unsigned int PBC     = 0x0;
const unsigned int APBC    = 0x1 << 3;
const unsigned int DIRC    = 0x1 << 4;
const unsigned int BC_MASK = PBC | APBC | DIRC;

const double INFTY = 10.0;

const int MAXDIMS = 4;