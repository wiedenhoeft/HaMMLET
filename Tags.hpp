#ifndef TAGS_HPP
#define TAGS_HPP


enum MappingType {combinations, independent};


// these empty classes are used as tags for template specialization, allowing the compiler to inline and optimize for different use cases

////////// These classes don't do anything by themselves, they are used as tags for template specialization //////////

class Dummy {};	// a tag for when the class is omitted, e.g. a transition matrix when no transitions are used

class Normal {};
class NormalParam {};
using NormalInverseGamma = NormalParam;
class NormalInverseGammaParam {};


class NormalVector {};
class NormalParamVector {};
using NormalInverseGammaVector = NormalParamVector;
class NormalInverseGammaParamVector {};


class InverseGamma {};
class InverseGammaParam {};
class NormalGammaParam {};


// using NormalGamma = NormalParam;

class Categorical {};
class CategoricalParam {};
using Dirichlet = CategoricalParam;
class DirichletParam {};


class CategoricalVector {};	// transition counts are SufficientStatistics<CategoricalVector>
class CategoricalParamVector {};
using DirichletVector = CategoricalParamVector;
class DirichletParamVector {};


class Geometric {};
class Beta {};
class BetaParam {};



// sampling tags for state sequence
class ForwardBackward {};	// forward-backward sampling
class Mixture {}; // sampling of each block individually
class DirectGibbs {};	// sample direct Gibbs, i.e. including transitions into and out of the state

////////// tags for data structures //////////
class Vector {}; // plain data structure for uncompressed sampling
class WaveletTree {};
class BreakpointArray {};





//////////  forward declarations //////////

template < typename StateSequenceType,
         typename EmissionDataStructure,
         typename EmissionDistType, // e.g. Normal
         typename ThetaDistType,	// e.g. NormalInverseGamma
         typename TransitionDistType, // e.g. DirichletVector
         typename InitialDistType,	// e.g. Dirichlet
         typename ThetaParamType,	// e.g. NormalInverseGammaParam
         typename TransitionParamType,	// e.g. DirichletParam
         typename InitialParamType 	// e.g. DirichletParam
         >
class HMM;


template < typename Type >
class StateSequence;

template < typename DataStructure,  typename DistType >
class Emissions;

template < typename DistType >
class Transitions;

template < typename DistType >
class Theta;

template < typename DistType >
class Initial;



template <typename T>
class Observation;

template <typename T>
class SufficientStatistics;

template <typename T>
class Distribution;

template <typename T>
class Conjugate;

template <typename DistType>
class ThetaHyperParam;


class Mapping;



#endif
