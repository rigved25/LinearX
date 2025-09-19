// include/linearx/partition/config.hpp
#pragma once
#include <linearx/partition/linear_partition.hpp>
#include <linearx/turbofold/linear_turbofold.hpp>

#define LP_TEMPLATE_TYPES \
    X(LinearPartition)    \
    X(TurboPartition)

#define LP_MODES(MACRO)           \
    MACRO(Mode::BEST)             \
    MACRO(Mode::PARTITION_INSIDE) \
    MACRO(Mode::PARTITION_OUTSIDE)
