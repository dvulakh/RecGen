
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
const gene& bp_domain::operator[](int i) const { return i ? genes.second : genes.first; }

// Comparators -- for map
/// Use pair comparators
bool bp_domain::operator==(const bp_domain& ot) const { return this->genes == ot.genes; }
bool bp_domain::operator<(const bp_domain& ot) const { return this->genes < ot.genes; }

/************************** BP MESSAGES ****************************/

// Constructor
/// Construct a new message with all elements initialized to nullval
bp_message::bp_message(long double nullval, int domain_sz)
{
	this->nullval = nullval;
	this->domain_sz = domain_sz;
}

// Internal marginals mutator
void bp_message::update_marginal(const gene &g, const long double &delta)
{
	auto it = this->marginals.find(g);
	if (it == this->marginals.end()) {
		it = this->marginals.insert({g, nullval * domain_sz}).first;
		this->unmapped.insert({g, domain_sz});
	}
	it->second += delta;
}

// Set & access mapped values
long double bp_message::operator[](const bp_domain& value) const
{
	auto it = this->probabilities.find(value);
	return it == this->probabilities.end() ? nullval : it->second;
}
long double bp_message::get_marginal(const gene &value)
{
	auto it = this->marginals.find(value);
	return it == this->marginals.end() ? nullval * domain_sz : it->second;
}
void bp_message::inc(const bp_domain& value, const long double delta)
{
	bool new_value = false;
	auto it = this->probabilities.find(value);
	if (it == this->probabilities.end()) {
		it = this->probabilities.insert({value, nullval}).first;
		new_value = true;
	}
	update_marginal(value[0], delta);
	update_marginal(value[1], delta);
	it->second += delta;
	this->unmapped.find(value[0])->second -= new_value;
	this->unmapped.find(value[1])->second -= new_value;
}
void bp_message::set(const bp_domain &value, const long double prob)
{ this->inc(value, prob - (*this)[value]); }

// Arithmetic operations
/// Add a message to this one
bp_message& bp_message::operator+=(bp_message& other)
{
	for (const auto &value : this->probabilities)
		this->inc(value.first, other[value.first]);
	for (const auto &value : other.probabilities)
		if (this->probabilities.find(value.first) == this->probabilities.end())
			this->inc(value.first, value.second);
	this->nullval += other.nullval;
	for (auto &value : this->marginals)
		value.second += other.nullval * this->unmapped.find(value.first)->second;
	return *this;
}
/// Multiply this message by another
bp_message& bp_message::operator*=(bp_message& other)
{
	for (const auto &value : this->probabilities)
		this->inc(value.first, (other[value.first] - 1) * value.second);
	for (const auto value : other.probabilities)
		if (this->probabilities.find(value.first) == this->probabilities.end())
			this->inc(value.first, (value.second - 1) * nullval);
	this->nullval *= other.nullval;
	return *this;
}
/// Multiply this message by a constant
bp_message& bp_message::operator*=(long double scalar)
{
	for (auto &value : this->probabilities)
		value.second *= scalar;
	for (auto &value : this->marginals)
		value.second *= scalar;
	this->nullval *= scalar;
	return *this;
}
bp_message& bp_message::operator/=(long double scalar) { return *this *= 1 / scalar; }

// Normalize this message for future use
void bp_message::normalize()
{
	/// Compute sum of probabilities
	long double sum = 0;
	for (auto value : this->probabilities)
		sum += value.second;
	sum += this->nullval * (domain_sz * domain_sz - this->probabilities.size());
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
