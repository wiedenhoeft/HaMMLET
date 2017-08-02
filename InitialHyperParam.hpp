#ifndef INITIALHYPERPARAM_HPP
#define INITIALHYPERPARAM_HPP

#include "includes.hpp"
#include "Tags.hpp"
#include "Conjugate.hpp"


// for the time being, tau is an alias for conjugates
template <typename ParamType>
using InitialHyperParam = Conjugate<ParamType>;


#endif
