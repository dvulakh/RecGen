
/********************************************************************
* Implements the interface provided by a struct representing the
* probability distribution estimates passed as BP 'messages'.
* Overloads arithmetic operators for ease of use.
********************************************************************/

#include "bp_message.h"

#include <algorithm>

/************************ DOMAIN ELEMENTS **************************/

// Constructor
/// Construct a domain element containing two genes
bp_domain::bp_domain(gene g1, gene g2)
{
	/// Use canonic sorted order
	this->genes = { std::min(g1, g2), std::max(g1, g2) };
}

// Accessor
/// Return the ith gene of the domain element
const gene& bp_domain::operator[](int i) { return i ? genes.second : genes.first; }

// Comparators -- for map
/// Use pair comparators
bool bp_domain::operator==(const bp_domain& ot) const { return this->genes == ot.genes; }
bool bp_domain::operator<(const bp_domain& ot) const { return this->genes < ot.genes; }

/************************** BP MESSAGES ****************************/

// Constructor
/// Construct a new message with all elements initialized to nullval
bp_message::bp_message(long double nullval) { this->nullval = nullval; }

// Set & access mapped values
long double& bp_message::operator[](const bp_domain& value)
{
	auto it = this->probabilities.find(value);
	if (it == this->probabilities.end())
		return this->probabilities.insert({ value, nullval }).first->second;
	return it->second;
}

// Arithmetic operations
/// Add a message to this one
bp_message& bp_message::operator+=(bp_message& other)
{
	for (auto &value : this->probabilities)
		value.second += other[value.first];
	for (auto value : other.probabilities)
		if (this->probabilities.find(value.first) == this->probabilities.end())
			(*this)[value.first] += value.second;
	this->nullval += other.nullval;
	return *this;
}
/// Multiply this message by another
bp_message& bp_message::operator*=(bp_message& other)
{
	for (auto &value : this->probabilities)
		value.second *= other[value.first];
	for (auto value : other.probabilities)
		if (this->probabilities.find(value.first) == this->probabilities.end())
			(*this)[value.first] *= value.second;
	this->nullval *= other.nullval;
	return *this;
}
/// Multiply this message by a constant
bp_message& bp_message::operator*=(long double scalar)
{
	for (auto &value : this->probabilities)
		value.second *= scalar;
	this->nullval *= scalar;
	return *this;
}
bp_message& bp_message::operator/=(long double scalar) { return *this *= 1 / scalar; }

// Normalize this message for future use
void bp_message::normalize(int domain_size)
{
	/// Compute sum of probabilities
	long double sum = 0;
	for (auto value : this->probabilities)
		sum += value.second;
	sum += this->nullval * (domain_size - this->probabilities.size());
	/// Divide by total sum
	*this /= sum;
}

// Extract maximum-probability domain element
bp_domain bp_message::extract_max()
{
	bp_domain max(1, 2);
	long double max_prop = 0;
	for (auto weight : this->probabilities)
		if (weight.second > max_prop)
			max = weight.first, max_prop = weight.second;
	return max;
}
