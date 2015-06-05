#ifndef DIGITH
#define DIGITH

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cctype>
#include <utility>
#include <algorithm>
#include <functional>

/*****************************************************************************/

namespace infinite_integer {

	class Digit;
	bool operator==(const Digit &_d1, const Digit &_d2);
	bool operator!=(const Digit &_d1, const Digit &_d2);
	bool operator<(const Digit &_d1, const Digit &_d2);
	bool operator>(const Digit &_d1, const Digit &_d2);
	bool operator<=(const Digit &_d1, const Digit &_d2);
	bool operator>=(const Digit &_d1, const Digit &_d2);
	Digit operator+(const Digit &_d1, const Digit &_d2);
	Digit operator-(const Digit &_d1, const Digit &_d2);
	Digit operator*(const Digit &_d1, const Digit &_d2);
	Digit operator/(const Digit &_d1, const Digit &_d2);
	Digit operator%(const Digit &_d1, const Digit &_d2);
	void swap(Digit &_d1, Digit &_d2);

	class Digit {
	private:
		using ll_t = long long;
		using ull_t = unsigned long long;

		friend bool operator==(const Digit &_d1, const Digit &_d2);
		friend bool operator<(const Digit &_d1, const Digit &_d2);
		friend void swap(Digit &_d1, Digit &_d2);

	public:
		// конструкторы
		Digit(ll_t _v = 0) : Digit(std::to_string(_v)) {}
		//Digit(const char *_v) : Digit(std::string(_v)) {}
		Digit(const std::string &_str);
		Digit(const Digit &_digit) = default;
		Digit(Digit &&_digit) throw();

		// операторы присвоения
		Digit &operator=(const Digit &_digit) = default;
		Digit &operator=(Digit &&_digit) throw();

		// преобразование в std::string
		operator std::string() const;

		// преобразование в long long (старший разряды обрезаются)
		// остаются 18 цифр
		long long to_long_long() const;

		// преобразование в bool;
		// false, если число равно нулю, независимо от знака;
		// true в любом другом случае
		explicit operator bool() const;

		// числовые операторы
		Digit &operator++();
		Digit &operator--();
		Digit operator++(int);
		Digit operator--(int);
		Digit &operator+=(const Digit &_digit);
		Digit &operator-=(const Digit &_digit);
		Digit &operator*=(const Digit &_digit);
		Digit &operator/=(const Digit &_digit);
		Digit &operator%=(const Digit &_digit);

		// обнулить
		Digit &set_to_zero();

		// изменить знак
		Digit &negate();

		// получить число в разряде _n
		unsigned get_digit_in_position(unsigned _n) const;
		// изменить число в разряде _n на _v
		void set_digit_in_position(unsigned _n, unsigned _v);
		// позиция самого большого ненулевого разряда (ln10)
		unsigned greatest_position() const;

		// модуль этого числа
		Digit abs() const;

		// удалить из вектора value неиспользуемые элементы (если старшие разряды = 0)
		void shrink_to_fit();

		// сдвиг всех цифр в числе на _u позиций влево (с добавлением нулей справа)
		Digit &shift_left(unsigned _n);
		// сдвиг всех цифр в числе на _u позиций вправо (с добавлением потерей символов справа)
		Digit &shift_right(unsigned _n);

		// нулевое значение
		static const Digit DIGIT_ZERO;

	private:
		// кол-во разрядов в ull_t
		static const unsigned SCALE = 18;	// 18
		// максимальное значение, которое хранит элемент массива value
		static const ull_t VALUE_ELEM_MAX = 999999999999999999;	// 999999999999999999
		// максимальное значение, которое может хранить ull_t
		static const ull_t ULL_T_MAX = ULLONG_MAX;

		// поэлементно сложить или вычесть векторы value (учитывая переполнение при сложении)
		void value_plus(const Digit &_digit);
		void value_minus(const Digit &_digit);	// (из большего меньшее)

		// умножить вектор value переданного числа на _u, учитывая переполение;
		// _u может быть равным 1 - 9;
		// если _u > 9, то производится умножение на (_u % 10);
		// для реалицации Digit::operator*=(const Digit &_digit)
		static Digit value_multi(const Digit &_digit, unsigned _u);

		// 10 в степени _u
		static Digit::ll_t pwr_10(unsigned _u);

		// сделать вектор value достаточно большим, чтобы его можно было складывать с _digit
		// (или находить разность)
		// для операторов += -=
		void reserve_enough_sum(const Digit &_digit);

		// сделать вектор value достаточно большим, чтобы его можно было умножать на _digit
		// для оператора *=
		void reserve_enough_multi(const Digit &_digit);

		// сделать вектор value достаточно большим,
		// чтобы его цифры можно было сдвинуть на _n разрядов влево
		void reserve_enough_shift_left(unsigned _n);

		// сделать вектор value достаточно большим, чтобы он мог содержать в позиции _n
		void reserve_enough_position(unsigned _n);

		// ? число = 0 (знак не имеет значения)
		bool is_zero() const;

		// само число. старшие разряды в элементах с большим индексом
		std::vector<ull_t> value;
		// знак числа
		bool positive;
	};

	/*****************************************************************************/

	const Digit Digit::DIGIT_ZERO(0);

	/*****************************************************************************/

	inline Digit::Digit(const std::string &_str) {
		positive = (_str.front() != '-');
		std::string str = _str.substr(_str.front() == '-' || _str.front() == '+' ? 1 : 0);
		// проверка корректности символов строки,
		// если найден неверный символ, строка образается до этого символа
		for (std::string::const_iterator b = str.cbegin(), e = str.cend(); b != e; ++b) {
			if (!std::isdigit(*b)) {
				str.erase(b, e);
				break;
			}
		}
		while (str.size()) {
			std::string temp;
			if (str.size() >= SCALE) {
				temp = str.substr(str.size() - SCALE);
				str.resize(str.size() - SCALE);
			}
			else {
				temp = std::move(str);
				str = "";
			}
			std::stringstream ss(temp);
			ull_t one_value;
			// преобразование строки в число
			ss >> one_value;
			value.emplace_back(one_value);
		}
	}

	// конструктор перемещения
	inline Digit::Digit(Digit &&_digit) throw()
		: value(std::move(_digit.value)), positive(_digit.positive) {}

	// оператор присваивания при перемещении
	inline Digit &Digit::operator=(Digit &&_digit) throw() {
		if (this != &_digit) {
			value = std::move(_digit.value);
			positive = _digit.positive;
		}
		return *this;
	}

	inline Digit &Digit::set_to_zero() {
		return *this = DIGIT_ZERO;
	}

	inline Digit &Digit::negate() {
		positive = (positive ? false : true);
		return *this;
	}

	inline Digit::operator std::string() const {
		if (this->is_zero()) {
			return "0";
		}
		std::ostringstream oss;
		if (!positive) {
			oss << '-';
		}
		std::vector<ull_t>::const_reverse_iterator iter = value.crend() - (greatest_position() / SCALE + 1);
		auto e = value.crend();
		oss << *iter++;
		// записываем элементы из value в oss, начиная со старших разрядов числа
		while (iter != e) {
			if (*iter == 0) {
				oss << std::string(SCALE, '0');
			}
			else {
				// кол-во нулей перед числом (элементом вектора value)
				unsigned zeros = 0;
				ull_t q = (VALUE_ELEM_MAX + 1) / 10;
				while (!(*iter / q)) {
					++zeros;
					q /= 10;
				}
				oss << std::string(zeros, '0');
				oss << *iter;
			}
			++iter;
		}
		return oss.str();
	}

	inline long long Digit::to_long_long() const {
		long long result = 0;
		for (unsigned u = 0; u != 18; ++u) {
			result += get_digit_in_position(u) * pwr_10(u);
		}
		if (!positive) {
			result *= -1;
		}
		return result;
	}

	inline Digit::operator bool() const {
		return !this->is_zero();
	}

	inline Digit &Digit::operator++() {
		return *this += 1;
	}

	inline Digit &Digit::operator--() {
		return *this -= 1;
	}

	inline Digit Digit::operator++(int) {
		Digit temp(*this);
		++*this;
		return temp;
	}

	inline Digit Digit::operator--(int) {
		Digit temp(*this);
		--*this;
		return temp;
	}

	inline Digit &Digit::operator+=(const Digit &_digit) {
		if (_digit.is_zero()) {
			return *this;
		}
		// чтобы все разряды точно влезли
		reserve_enough_sum(_digit);
		if (is_zero()) {
			return *this = _digit;
		}
		if (positive == _digit.positive) {
			// если оба слагаемых имеют один знак, то складываем их модули,
			// а знак оставляем тот же
			value_plus(_digit);
		}
		else {
			// если знаки разные, то определяем знак результата
			// и вычитаем из модуля большего числа модуль меньшего
			bool result_positive = (this->abs() > _digit.abs()) == positive;
			Digit min_summand = std::min(this->abs(), _digit.abs());
			Digit max_summand = std::max(this->abs(), _digit.abs());
			max_summand.value_minus(min_summand);
			*this = std::move(max_summand);
			positive = result_positive;
		}
		return *this;
	}

	inline Digit &Digit::operator-=(const Digit &_digit) {
		return *this += Digit(_digit).negate();
	}

	inline Digit &Digit::operator*=(const Digit &_digit) {
		// проверка умножения на ноль
		if (is_zero() || _digit.is_zero()) {
			return *this = DIGIT_ZERO;
		}
		// знак результата
		bool result_positive = (_digit.positive == positive);
		// проверка умножения на 1 и -1
		if (_digit.abs() == 1) {
			positive = result_positive;
			return *this;
		}
		// а, может, this это 1 или -1...
		if (this->abs() == 1) {
			*this = _digit;
			positive = result_positive;
			return *this;
		}
		// чтобы все разряды точно влезли
		reserve_enough_multi(_digit);
		// результат
		Digit temp = DIGIT_ZERO;
		// номер наибольшего разряда в _digit
		int digit_pos = _digit.greatest_position();
		// умножаем this на цифру из каждой позиции из _digit с соответствующий сдвигом
		while (digit_pos >= 0) {
			temp += value_multi(*this, _digit.get_digit_in_position(digit_pos)).shift_left(digit_pos);
			--digit_pos;
		}
		*this = std::move(temp);
		positive = result_positive;
		return *this;
	}

	inline Digit &Digit::operator/=(const Digit &_digit) {
		if (_digit.is_zero()) {
			throw std::runtime_error("! division by zero");
		}
		// знак результата
		bool result_positive = (_digit.positive == positive);
		// проверка деления на 1 и -1
		if (_digit.abs() == 1) {
			positive = result_positive;
			return *this;
		}
		// а, может, this меньше _digit и результат просто ноль...
		if (this->abs() < _digit.abs()) {
			return *this = DIGIT_ZERO;
		}
		// или делимое и делитель равны по модулю...
		if (this->abs() == _digit.abs()) {
			*this = 1;
			positive = result_positive;
			return *this;
		}
		Digit d1 = this->abs();					// делимое
		Digit d2 = _digit.abs();				// делитель
		Digit temp = DIGIT_ZERO;				// частное (результат)
		// имитация деления столбиком (делим только модули, знак уже известен)
		int shift = this->greatest_position() - _digit.greatest_position();
		while (d1 >= d2) {
			d2.shift_left(shift);
			while (d2 > d1) {
				d2.shift_right(1);
				--shift;
			}
			Digit d2_temp = d2;
			int one_digit_temp = 1;
			while (d1 >= d2) {
				d2 += d2_temp;
				++one_digit_temp;
			}
			d2 -= d2_temp;
			--one_digit_temp;
			d1 -= d2;
			d2 = _digit.abs();
			temp.set_digit_in_position(shift, one_digit_temp);
			--shift;
		}
		*this = std::move(temp);
		positive = result_positive;
		return *this;
	}

	inline Digit &Digit::operator%=(const Digit &_digit) {
		if (_digit.is_zero()) {
			throw std::runtime_error("! division by zero");
		}
		if (_digit == 1) {
			return *this = DIGIT_ZERO;
		}
		if (_digit == 2) {
			return *this = (get_digit_in_position(0) % 2);
		}
		return *this -= (*this / _digit) * _digit;
	}

	inline unsigned Digit::get_digit_in_position(unsigned _n) const {
		if (_n >= value.size() * SCALE) {
			return 0;
		}
		std::vector<ull_t>::size_type i = _n / SCALE;
		std::vector<ull_t>::size_type j = _n % SCALE;
		ull_t one_value = *(value.cbegin() + i);
		one_value /= pwr_10(j);
		one_value %= 10;
		return static_cast<int>(one_value);
	}

	inline void Digit::set_digit_in_position(unsigned _n, unsigned _v) {
		// если требуемого разряда еще нет, он создается
		reserve_enough_position(_n);
		std::vector<ull_t>::size_type i = _n / SCALE;
		std::vector<ull_t>::size_type j = _n % SCALE;
		// элемент вектора, в котором находится изменяемый разряд
		ull_t &one_value = *(value.begin() + i);
		// то, что раньше было в этом разряде
		int prev = (one_value / pwr_10(j)) % 10;
		// разница между новым и старым значением
		int diff = _v - prev;
		// внесение изменений
		one_value += diff * pwr_10(j);
	}

	inline unsigned Digit::greatest_position() const {
		int ind = value.size() * SCALE;
		while (get_digit_in_position(--ind) == 0);
		return ind;
	}

	inline Digit Digit::abs() const {
		Digit temp(*this);
		temp.positive = true;
		return temp;
	}

	inline void Digit::shrink_to_fit() {
		std::vector<ull_t>::size_type ind = value.size() - 1;
		// удалить элементы, содержащие старшие разряды, если они равны 0
		while (ind != 0 && value[ind] == 0) {
			--ind;
		}
		value.resize(ind + 1);
		value.shrink_to_fit();
	}

	inline void Digit::value_plus(const Digit &_digit) {
		for (std::vector<ull_t>::size_type ind = 0, max_ind = value.size(),
			max_ind_digit = _digit.value.size(); ind != max_ind_digit; ++ind) {
			// складываем числа поэлементно в векторах value
			if ((value[ind] += _digit.value[ind]) > VALUE_ELEM_MAX) {
				// произошло переполнение элемента вектора value
				value[ind] -= VALUE_ELEM_MAX + 1;
				if (ind != max_ind - 1) {
					// если элемент не самый старший, просто увеличим более старший на 1
					++value[ind + 1];
				}
				else {
					// если элемент самый старший, добавим новый элемент со значением 1
					value.emplace_back(1);
				}
			}
		}
	}

	inline void Digit::value_minus(const Digit &_digit) {
		for (std::vector<ull_t>::size_type ind = 0, max_ind = value.size(),
			max_ind_digit = _digit.value.size(); ind != max_ind_digit; ++ind) {
			// находим поэлементную разность элементов вектора value
			if ((value[ind] -= _digit.value[ind]) > VALUE_ELEM_MAX) {
				// произошло переполнение элемента вектора value
				value[ind] *= -1;
				value[ind] = VALUE_ELEM_MAX - value[ind] + 1;
				if (ind != max_ind - 1) {
					// если элемент не самый старший, просто уменьшим более старший на 1
					--value[ind + 1];
				}
				else {
					// мы вычитаем из большего меньшее,
					// поэтому самый старший разряд не может переполниться
					;
				}
			}
		}
	}

	inline Digit Digit::value_multi(const Digit &_digit, unsigned _u) {
		if (_u > 9) {
			// чтобы корректно работало
			_u %= 10;
		}
		if (_u == 0) {
			return DIGIT_ZERO;
		}
		Digit temp = _digit;
		// то, на сколько нужно будет увеличить следующий элемент вектора value при переполнении
		ull_t next_val = 0;
		for (std::vector<Digit::ull_t>::size_type ind = 0,
			max_ind = temp.value.size(); ind != max_ind; ++ind) {
			// поэлементную умножаем вектор value на _u
			temp.value[ind] *= _u;
			temp.value[ind] += next_val;
			next_val = 0;
			if (temp.value[ind] > VALUE_ELEM_MAX) {
				// произошло переполнение элемента вектора value
				next_val = temp.value[ind] / (VALUE_ELEM_MAX + 1);
				// вернем значение текущего элемента вектора в нормальное состояние
				temp.value[ind] %= VALUE_ELEM_MAX + 1;
				if (ind != max_ind - 1) {
					; // если элемент не самый старший, то он будет увеличен во время след. итерации
				}
				else {
					// если элемент самый старший, добавим новый элемент со нужным значением
					temp.value.emplace_back(next_val);
				}
			}
		}
		return temp;
	}

	inline Digit::ll_t Digit::pwr_10(unsigned _u) {
		Digit::ll_t result = 1;
		while (_u--) {
			result *= 10;
		}
		return result;
	}

	inline Digit &Digit::shift_left(unsigned _n) {
		// ноль сдвигать бессмысленно, и на 0 позиций тоже
		if (this->is_zero() || _n == 0) {
			return *this;
		}
		// найдем наибольший ненулевой разряд числа
		int ind = greatest_position();
		// сделаем value достаточно большим
		reserve_enough_shift_left(_n);
		// выполняем сдвиг разрядов, включая нулевой
		while (ind >= 0) {
			set_digit_in_position(ind + _n, get_digit_in_position(ind));
			--ind;
		}
		// заполняем нулями разряды справа
		for (ind = 0; ind != _n; ++ind) {
			set_digit_in_position(ind, 0);
		}
		return *this;
	}

	inline Digit &Digit::shift_right(unsigned _n) {
		if (this->is_zero() || _n == 0) {
			return *this;
		}
		// максимальный ненулевой разряд числа до сдвига
		unsigned greatest_position_before = greatest_position();
		if (_n > greatest_position_before) {
			// если сдвиг на больше позиций, чем есть в числе
			return *this = DIGIT_ZERO;
		}
		// выполняем сдвиг
		for (unsigned ind = 0; ind <= greatest_position_before; ++ind) {
			set_digit_in_position(ind, get_digit_in_position(ind + _n));
		}
		return *this;
	}

	inline void Digit::reserve_enough_sum(const Digit &_digit) {
		if (value.size() < _digit.value.size()) {
			value.resize(_digit.value.size());
		}
	}

	inline void Digit::reserve_enough_multi(const Digit &_digit) {
		value.resize((greatest_position() + _digit.greatest_position()) / SCALE + 1);
	}

	inline void Digit::reserve_enough_shift_left(unsigned _n) {
		reserve_enough_position(greatest_position() + _n);
	}

	inline void Digit::reserve_enough_position(unsigned _n) {
		if (_n >= value.size() * SCALE) {
			value.resize((_n / SCALE) + 1);
		}
	}

	inline bool Digit::is_zero() const {
		return std::find_if(value.cbegin(), value.cend(),
			std::bind(std::not_equal_to<ull_t>(), std::placeholders::_1, 0)) == value.cend();
	}

	/*****************************************************************************/

	// IO

	inline std::ostream &operator<<(std::ostream &_out, const Digit &_digit) {
		return _out << std::string(_digit);
	}

	inline std::istream &operator>>(std::istream &_in, Digit &_digit) {
		std::string temp;
		_in >> temp;
		_digit = temp;
		return _in;
	}

	// операторы с 2 симметричными параметрами

	inline bool operator==(const Digit &_d1, const Digit &_d2) {
		// положительный и отрицательный нули равны
		if (_d1.is_zero() && _d2.is_zero()) {
			return true;
		}
		// числа с разными знаками не равны
		if (_d1.positive != _d2.positive) {
			return false;
		}
		auto not_zero = std::bind(std::not_equal_to<Digit::ull_t>(), std::placeholders::_1, 0);
		// разность длины векторов value этого числа (this) и _digit
		std::vector<int>::difference_type size_diff = _d1.value.size() - _d2.value.size();
		// если вектора имеют размеры и у большего есть числа старшие разряды не равные 0, то числа не равны
		if (size_diff > 0) {
			if (std::find_if(_d1.value.crbegin(), _d1.value.crbegin() + size_diff, not_zero) !=
				_d1.value.crbegin() + size_diff) {
				return false;
			}
		}
		if (size_diff < 0) {
			if (std::find_if(_d2.value.crbegin(), _d2.value.crbegin() - size_diff, not_zero) !=
				_d2.value.crbegin() - size_diff) {
				return false;
			}
		}
		// теперь сравним те разряды, которые есть у обоих чисел, учитывая знаки чисел
		return _d1.positive == _d2.positive &&
			std::equal(_d1.value.cbegin(),
			(_d1.value.cbegin() + std::min(_d1.value.size(), _d2.value.size())), _d2.value.cbegin());
	}


	inline bool operator!=(const Digit &_d1, const Digit &_d2) {
		return !(_d1 == _d2);
	}

	inline bool operator<(const Digit &_d1, const Digit &_d2) {
		// положительный и отрицательный нули равны, и ни один из них не меньше другого
		if (_d1.is_zero() && _d2.is_zero()) {
			return false;
		}
		if (_d1.positive ^ _d2.positive) {
			return !_d1.positive;
		}
		auto not_zero = std::bind(std::not_equal_to<Digit::ull_t>(), std::placeholders::_1, 0);
		// разность длины векторов value этого числа (this) и _digit
		std::vector<int>::difference_type size_diff = _d1.value.size() - _d2.value.size();
		// если вектора имеют размеры и у большего есть числа старшие разряды не равные 0, то одно число меньше другого
		if (size_diff > 0) {
			if (std::find_if(_d1.value.crbegin(), _d1.value.crbegin() + size_diff, not_zero) !=
				_d1.value.crbegin() + size_diff) {
				return !_d1.positive;
			}
		}
		if (size_diff < 0) {
			if (std::find_if(_d2.value.crbegin(), _d2.value.crbegin() - size_diff, not_zero) !=
				_d2.value.crbegin() - size_diff) {
				return _d1.positive;
			}
		}
		// минимум среди размеров сравниваемых вектров
		std::vector<int>::difference_type common_size = std::min(_d1.value.size(), _d2.value.size());
		// теперь сравним те разряды, которые есть у обоих чисел, учитывая знаки чисел
		return std::lexicographical_compare(
			_d1.value.crend() - common_size, _d1.value.crend(),
			_d2.value.crend() - common_size, _d2.value.crend())
			== _d1.positive;
	}


	inline bool operator>(const Digit &_d1, const Digit &_d2) {
		return !(_d1 < _d2 || _d1 == _d2);
	}

	inline bool operator<=(const Digit &_d1, const Digit &_d2) {
		return !(_d1 > _d2);
	}

	inline bool operator>=(const Digit &_d1, const Digit &_d2) {
		return !(_d1 < _d2);
	}

	inline Digit operator+(const Digit &_d1, const Digit &_d2) {
		return Digit(_d1) += _d2;
	}

	inline Digit operator-(const Digit &_d1, const Digit &_d2) {
		return Digit(_d1) -= _d2;
	}

	inline Digit operator*(const Digit &_d1, const Digit &_d2) {
		return Digit(_d1) *= _d2;
	}

	inline Digit operator/(const Digit &_d1, const Digit &_d2) {
		return Digit(_d1) /= _d2;
	}

	inline Digit operator%(const Digit &_d1, const Digit &_d2) {
		return Digit(_d1) %= _d2;
	}

	inline void swap(Digit &_d1, Digit &_d2) throw() {
		using std::swap;
		swap(_d1.value, _d2.value);
		swap(_d1.positive, _d2.positive);
	}

}	// namespace infinite_integer end

#endif
