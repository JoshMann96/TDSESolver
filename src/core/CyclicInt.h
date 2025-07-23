/**
 * @file CyclicInt.h
 * @brief Cyclic integer class for managing cyclic indices. Useful for managing a finite history of values without moving data around to keep things in order.
 */
#pragma once

/// A cyclic integer class that wraps around a maximum value. Useful for managing the local wavefunction and potential history.
template<typename T>
class cyclic_int
{
	static_assert(std::is_integral<T>::value, "cyclic_int only works with integral types");
protected:
	T val, max;
public:
    /// Default constructor initializes the value to 0 and the maximum to 0.
	cyclic_int() : val(0), max(0) {};

	/**
	 * Constructor initializes the value to 0 and the maximum to the given value.
	 * @param max The maximum value (exclusive) for the cyclic integer.
	 */
	cyclic_int(T max) : val(0), max(max) {};

	/**
	 * Constructor initializes the value to the given value and the maximum to the given value.
	 * @param val The initial value of the cyclic integer.
	 * @param max The maximum value (exclusive) for the cyclic integer.
	 */
	cyclic_int(T val, T max) : val(val), max(max) {};

	inline void increment() { val = (val + 1) % max; };
	inline cyclic_int& operator++() { increment(); return *this; }; //prefix
	inline cyclic_int operator++(int) { cyclic_int c = *this; increment(); return c; }; //postfix
	inline cyclic_int operator+(T n) { cyclic_int c(max); c.val = (val + n) % max; return c; };
	inline cyclic_int operator+(int n) { cyclic_int c(max); c.val = (val + n) % max; return c; };
	inline cyclic_int& operator+=(T n) { val = (val + n) % max; return *this; };
	inline cyclic_int& operator=(T n) { val = n % max; return *this; };
	inline cyclic_int operator-(T n) { 
		if(n >= max) throw std::out_of_range("cyclic_int: Cannot subtract more than max value");
		cyclic_int c(max); c.val = (val + max - n) % max; return c;
	};
	inline cyclic_int operator-(int n) {
		if(n >= max) throw std::out_of_range("cyclic_int: Cannot subtract more than max value");
		cyclic_int c(max); c.val = (val + max - n) % max; return c;
	};
	inline operator size_t() const { return val; };
};