#include "InteractionPlatelet.hpp"

Interaction::Interaction() : i(0), j(0), fn(0.0), ft(0.0) {}
Interaction::Interaction(size_t I, size_t J, size_t KI, size_t KJ)
    : i(I), j(J), ki(KI), kj(KJ), fn(0.0), ft(0.0), delta_t(0.0) {}
