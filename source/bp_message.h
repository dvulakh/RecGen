
/********************************************************************
* Defines the interface provided by a struct representing the
* probability distribution estimates passed as BP 'messages'.
* Overloads arithmetic operators for ease of use.
********************************************************************/

#ifndef BP_MESSAGE_H
#define BP_MESSAGE_H

#include "poisson_pedigree.h"

#include <map>

// Use a canonic form of gene pairs as domain elements of the marginals
// for BP messages
struct bp_domain
{
private:
	std::pair<gene, gene> genes;
public:
	// Constructor
	/// Construct a domain element containing two genes
	bp_domain(gene g1, gene g2);
	// Accessor
	/// Return the ith gene of the domain element
	const gene& operator[](int i) const;
	// Comparators -- for map
	bool operator==(const bp_domain& ot) const;
	bool operator<(const bp_domain& ot) const;
};

// Representation of messages for BP algorithm
struct bp_message
{
private:
	// Private members
	/// Represent non-default values with a map
	std::map<bp_domain, long double> probabilities;
	/// Represent total probabilities per half-value
	std::map<gene, long double> marginals;
	/// Represent number of unassigned domain elements
	std::map<gene, int> unmapped;
	/// Set default value at initialization time
	long double nullval;
	/// Set domain size (per gene) at initialization time
	int domain_sz;
	// Internal mutation functions
	void update_marginal(const gene& g, const long double& delta);
public:
	// Constructor
	/// Construct a new message with all elements initialized to nullval
	bp_message(long double nullval, int domain_sz);
	// Set & access mapped values
	long double operator[](const bp_domain& value) const;
	long double get_marginal(const gene& value);
	void inc(const bp_domain& value, const long double delta);
	void set(const bp_domain& value, const long double prob);
	// Arithmetic operations
	/// Add a message to this one
	bp_message& operator+=(bp_message& other);
	/// Multiply this message by another
	bp_message& operator*=(bp_message& other);
	/// Multiply this message by a constant
	bp_message& operator*=(long double scalar);
	/// Divide this message by a constant
	bp_message& operator/=(long double scalar);
	// Normalize this message for future use
	bp_message& normalize();
	// Purge pairs from memory, keeping only marginals
	bp_message& purge();
	// Extract maximum value
	bp_domain extract_max();
};

#endif
