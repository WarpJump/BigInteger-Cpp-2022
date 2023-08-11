#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using std::string;

class BigInteger {
  static const long long kBase = 1e9;
  static const int kBaseModule = 9;

  std::vector<long long> polynom_digits;
  unsigned int size_;
  bool isNonNegative_;

  void deleteInsignificantZeros();

  void kBaseFormRestoration();

  string addInsignificantZeros(const string& incomplete) const {
    int numberOfZeros = kBaseModule - static_cast<int>(incomplete.length());
    if (numberOfZeros > 0) {
      return string(static_cast<size_t>(numberOfZeros), '0') + incomplete;
    }
    return incomplete;
  }

  BigInteger collapseDegree(int) const;

 public:
  friend bool operator<(const BigInteger&, const BigInteger&);
  friend bool operator==(const BigInteger&, const BigInteger&);

  BigInteger() : polynom_digits(0), size_(1), isNonNegative_(true) {}

  BigInteger(const int constructor) : size_(1), isNonNegative_(constructor >= 0) {
    polynom_digits.push_back(abs(constructor));
  }

  explicit BigInteger(const string& str) : isNonNegative_(true) {
    int length = static_cast<int>(str.length());
    int key = 0;
    if (str[0] == '-') {
      isNonNegative_ = false;
      ++key;
    }
    size_ = (length - key - 1) / kBaseModule + 1;

    for (int BigInteger_digit = length - 1; BigInteger_digit >= key;
         BigInteger_digit -= kBaseModule) {
      long long basic_sub_int = 0;
      for (int runner = kBaseModule - 1; runner >= 0; --runner) {
        if (BigInteger_digit - runner >= key) {
          basic_sub_int *= 10;
          basic_sub_int += (str[BigInteger_digit - runner] - '0');
        }
      }
      polynom_digits.push_back(basic_sub_int);
    }
  }

  BigInteger restoreDegree(int);

  int size() const { return size_; }
  bool sign() const { return isNonNegative_; }

  BigInteger operator-() const {
    BigInteger negative = *this;
    negative.isNonNegative_ = !negative.isNonNegative_;
    return negative;
  }

  explicit operator bool() const {
    return !(size_ == 1 && polynom_digits[0] == 0);
  }

  string toString() const;

  BigInteger& operator+=(const BigInteger&);

  BigInteger& operator-=(const BigInteger&);

  BigInteger& operator*=(const BigInteger&);

  BigInteger& operator%=(long long divider) {
    *this = abs(polynom_digits[0] %= divider);
    isNonNegative_ = isNonNegative_ == (divider >= 0);
    return *this;
  }

  BigInteger& operator%=(const BigInteger&);

  BigInteger& operator/=(const BigInteger&);

  BigInteger operator--(int) {
    BigInteger temp = *this;
    *this -= 1;
    return temp;
  }

  BigInteger operator++(int) {
    BigInteger temp = *this;
    *this += 1;
    return temp;
  }

  BigInteger& operator--() {
    *this -= 1;
    return *this;
  }

  BigInteger& operator++() {
    *this += 1;
    return *this;
  }
};

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
  BigInteger ans = first;
  return ans *= second;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
  BigInteger ans = first;
  return ans += second;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
  BigInteger ans = first;
  return ans -= second;
}

BigInteger operator/(const BigInteger& first, const BigInteger& divider) {
  BigInteger answer = first;
  answer /= divider;
  return answer;
}

BigInteger operator%(const BigInteger& first, const BigInteger& divider) {
  BigInteger ans = first - (first / divider) * divider;
  return ans;
}

bool operator<(const BigInteger& first, const BigInteger& second) {
  if (first.isNonNegative_ != second.isNonNegative_) {
    return (first.polynom_digits[0] != 0) && (first.isNonNegative_ == false);
  }

  if (first.size_ != second.size_) {
    return (first.size_ < second.size_) == first.isNonNegative_;
  }

  for (int i = first.size_ - 1; i >= 0; --i) {
    if (first.polynom_digits[i] != second.polynom_digits[i]) {
      return (first.polynom_digits[i] < second.polynom_digits[i]) ==
             first.isNonNegative_;
    }
  }
  return false;
}

bool operator==(const BigInteger& first, const BigInteger& second) {
  if (first.size_ != second.size_) {
    return false;
  }
  if (first.isNonNegative_ != second.isNonNegative_) {
    return (first.polynom_digits[0] == 0 && second.polynom_digits[0] == 0);
  }

  for (int i = first.size_ - 1; i >= 0; --i) {
    if (first.polynom_digits[i] != second.polynom_digits[i]) return false;
  }
  return true;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
  return !(first == second);
}

bool operator>(const BigInteger& first, const BigInteger& second) {
  return second < first;
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
  return !(first > second);
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
  return !(first < second);
}

BigInteger BigInteger::collapseDegree(int degree) const {
  BigInteger collapsed;
  collapsed.size_ = size_ - degree;
  collapsed.polynom_digits.resize(collapsed.size_ + 1);
  if (size_ > static_cast<size_t>(degree)) {
    std::memcpy(&(collapsed.polynom_digits[0]), &(polynom_digits[degree]),
                sizeof(long long) * (size_ - degree));
  }
  collapsed.deleteInsignificantZeros();
  return collapsed;
}

BigInteger BigInteger::restoreDegree(int zeros) {
  BigInteger temp = *this;
  temp.polynom_digits.resize(size_ + zeros);

  std::memset(&(temp.polynom_digits[0]), 0, sizeof(long long) * zeros);
  std::memcpy(&(temp.polynom_digits[zeros]), &(polynom_digits[0]),
              sizeof(long long) * (size_));
  temp.size_ = size_ + zeros;
  *this = temp;
  return *this;
}

void BigInteger::deleteInsignificantZeros() {
  size_ = polynom_digits.size();
  if (size_ == 0) {
    polynom_digits.push_back(0);
    ++size_;
  }
  while (polynom_digits[size_ - 1] == 0 && size_ > 1) {
    polynom_digits.pop_back();
    --size_;
  }
}

void BigInteger::kBaseFormRestoration() {
  long long overdigit = 0;
  for (unsigned int i = 0; i < size_; ++i) {
    polynom_digits[i] += overdigit;
    overdigit = 0;
    if (polynom_digits[i] > kBase) {
      overdigit += polynom_digits[i] / kBase;
      polynom_digits[i] %= kBase;
    } else if (polynom_digits[i] < 0) {
      if (polynom_digits[i] <= -kBase) {
        overdigit = polynom_digits[i] / kBase;
        polynom_digits[i] = polynom_digits[i] % kBase;
      }
      if (polynom_digits[i] < 0) {
        --overdigit;
        polynom_digits[i] += kBase;
      }
    }
  }
  if (overdigit != 0) {
    polynom_digits.push_back(overdigit);
  }
  size_ = polynom_digits.size();
  deleteInsignificantZeros();
}

string BigInteger::toString() const {
  int size = polynom_digits.size();
  string str = std::to_string(polynom_digits[size - 1]);
  if ((isNonNegative_ == false) && bool(*this)) {
    str = '-' + str;
  }

  for (int i = size - 2; i >= 0; --i) {
    str += addInsignificantZeros(std::to_string(polynom_digits[i]));
  }
  return str;
}

BigInteger& BigInteger::operator+=(const BigInteger& summand) {
  if (summand.isNonNegative_ != isNonNegative_) {
    *this -= (-summand);
    return *this;
  }

  unsigned int runner = 0;
  while (runner < summand.size_) {
    if (runner < size_) {
      polynom_digits[runner] += summand.polynom_digits[runner];
    } else {
      polynom_digits.push_back(summand.polynom_digits[runner]);
    }
    ++runner;
  }
  kBaseFormRestoration();
  return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& subtrahend) {
  if (subtrahend.isNonNegative_ != isNonNegative_) {
    return *this += (-subtrahend);
  }
  if ((*this < subtrahend) == isNonNegative_) {
    isNonNegative_ = !isNonNegative_;
  }

  unsigned int runner = 0;
  while (runner < subtrahend.size_) {
    if (runner < size_) {
      polynom_digits[runner] -= subtrahend.polynom_digits[runner];
      if (isNonNegative_ != subtrahend.isNonNegative_) {
        polynom_digits[runner] *= -1;
      }
    } else {
      polynom_digits.push_back(subtrahend.polynom_digits[runner]);
    }
    ++runner;
  }

  kBaseFormRestoration();
  return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& multiplier) {
  BigInteger answer;
  answer.polynom_digits.resize(size_ + multiplier.size_);
  answer.size_ = size_ + multiplier.size_;
  answer.isNonNegative_ = (isNonNegative_ == multiplier.isNonNegative_);

  long long digit;
  for (unsigned int i = 0; i < size_; ++i) {
    for (unsigned int j = 0; j < multiplier.size_; ++j) {
      digit = polynom_digits[i] * multiplier.polynom_digits[j];
      answer.polynom_digits[i + j] += digit % kBase;
      answer.polynom_digits[i + j + 1] += digit / kBase;
    }
  }
  answer.kBaseFormRestoration();
  *this = answer;
  return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& divider) {
  *this -= (*this / divider) * divider;
  return *this;
}

BigInteger& BigInteger::operator/=(const BigInteger& divider) {
  if (size_ < divider.size_) {
    *this = BigInteger(0);
    return *this;
  }
  BigInteger answer;
  BigInteger approximation;
  answer.size_ = size_ - divider.size_ + 1;
  answer.polynom_digits.resize(answer.size_);
  answer.isNonNegative_ = (isNonNegative_ == divider.isNonNegative_);
  isNonNegative_ = true;

  for (int degree = answer.size_ - 1; degree >= 0; --degree) {
    BigInteger collapsed = collapseDegree(degree);
    collapsed.isNonNegative_ = true;
    long long lowerLimit;
    long long upperLimit;
    long long average;
    if (collapsed.size_ > divider.size_) {
      upperLimit = (collapsed.polynom_digits[collapsed.size_ - 1] * kBase +
                    collapsed.polynom_digits[collapsed.size_ - 2]) /
                       divider.polynom_digits[divider.size_ - 1] +
                   1;
    } else {
      upperLimit = collapsed.polynom_digits[collapsed.size_ - 1] /
                       divider.polynom_digits[divider.size_ - 1] +
                   1;
    }
    lowerLimit = upperLimit >> 2;

    while (lowerLimit < upperLimit - 1) {
      average = (lowerLimit + upperLimit) >> 1;
      approximation = divider * average;
      approximation.isNonNegative_ = true;

      if (approximation <= collapsed) {
        lowerLimit = average;
      } else {
        upperLimit = average;
      }
    }

    answer.polynom_digits[degree] = lowerLimit;

    approximation = divider * lowerLimit;
    approximation.restoreDegree(degree);
    approximation.isNonNegative_ = true;

    *this -= approximation;
  }
  answer.kBaseFormRestoration();
  *this = answer;
  return *this;
}

class Rational {
  static const int kBaseModule = 9;
  static const int kDoublePrecision = 308;

  BigInteger numerator_;
  BigInteger denominator_;
  bool isNonNegative_;

  Rational simplifyFraction();

  static BigInteger greatestCommonDivisor(BigInteger, BigInteger);

 public:
  friend bool operator<(const Rational&, const Rational&);
  friend bool operator==(const Rational&, const Rational&);
  Rational() : numerator_(0), denominator_(1), isNonNegative_(true) {}

  Rational(const BigInteger& integer)
      : denominator_(1), isNonNegative_(integer.sign()) {
    numerator_ = isNonNegative_ ? integer : -integer;
  }

  Rational(double value)
      : numerator_(abs(static_cast<int>(value))),
        denominator_(1),
        isNonNegative_(value >= 0) {}

  string toString() const;

  string asDecimal(size_t precision = 0) const;

  explicit operator double() const {
    double value = atof(asDecimal(kDoublePrecision).data());
    return value;
  }

  Rational operator-() const {
    Rational temp = *this;
    temp.isNonNegative_ = !(temp.isNonNegative_);
    return temp;
  }

  Rational& operator*=(const Rational& number) {
    numerator_ *= number.numerator_;
    denominator_ *= number.denominator_;
    isNonNegative_ = (isNonNegative_ == number.isNonNegative_);
    return *this;
  }

  Rational& operator/=(const Rational& number) {
    numerator_ *= number.denominator_;
    denominator_ *= number.numerator_;
    isNonNegative_ = (isNonNegative_ == number.isNonNegative_);
    return *this;
  }

  Rational& operator-=(const Rational& number) {
    if (isNonNegative_ != number.isNonNegative_) {
      return *this += (-number);
    }
    numerator_ *= number.denominator_;
    BigInteger subtrahend(number.numerator_);
    subtrahend *= denominator_;

    if (subtrahend > numerator_) {
      subtrahend -= numerator_;
      numerator_ = subtrahend;
      isNonNegative_ = !isNonNegative_;
    } else {
      numerator_ -= subtrahend;
    }
    denominator_ *= number.denominator_;
    simplifyFraction();
    return *this;
  }

  Rational& operator+=(const Rational& number) {
    if (isNonNegative_ != number.isNonNegative_) {
      return *this -= (-number);
    }
    numerator_ =
        numerator_ * number.denominator_ + number.numerator_ * denominator_;
    denominator_ *= number.denominator_;
    simplifyFraction();
    return *this;
  }
};

Rational operator*(const Rational& first, const Rational& second) {
  Rational ans = first;
  return ans *= second;
}

Rational operator/(const Rational& first, const Rational& second) {
  Rational ans = first;
  return ans /= second;
}

Rational operator+(const Rational& first, const Rational& second) {
  Rational ans = first;
  return ans += second;
}

Rational operator-(const Rational& first, const Rational& second) {
  Rational ans = first;
  return ans -= second;
}

bool operator<(const Rational& first, const Rational& second) {
  if (first.isNonNegative_ != second.isNonNegative_) {
    if (!bool(first.numerator_) && !bool(second.numerator_)) {
      return false;
    }
    return first.isNonNegative_ == false;
  }
  return ((first.numerator_ * second.denominator_) <
              (second.numerator_ * first.denominator_) ==
          first.isNonNegative_);
}

bool operator==(const Rational& first, const Rational& second) {
  if (first.isNonNegative_ != second.isNonNegative_) {
    return (!bool(first.numerator_) && !bool(second.numerator_));
  }
  return ((first.numerator_ * second.denominator_) ==
          (second.numerator_ * first.denominator_));
}

bool operator>(const Rational& first, const Rational& second) {
  return second < first;
}

bool operator<=(const Rational& first, const Rational& second) {
  return !(first > second);
}

bool operator>=(const Rational& first, const Rational& second) {
  return !(first < second);
}

bool operator!=(const Rational& first, const Rational& second) {
  return !(first == second);
}

BigInteger Rational::greatestCommonDivisor(BigInteger first,
                                           BigInteger second) {
  BigInteger temporary(second);
  while (bool(first) && bool(second)) {
    second = first % second;
    first = temporary;
    temporary = second;
  }
  if (!bool(first)) {
    return second;
  }
  return first;
}

Rational Rational::simplifyFraction() {
  BigInteger gcd = greatestCommonDivisor(numerator_, denominator_);
  numerator_ /= gcd;
  denominator_ /= gcd;
  return *this;
}

string Rational::toString() const {
  Rational number = *this;
  number.simplifyFraction();
  string stringRepresentation;
  if (!number.isNonNegative_) {
    stringRepresentation += '-';
  }
  stringRepresentation += number.numerator_.toString();
  if (number.denominator_ > BigInteger(1)) {
    stringRepresentation += "/" + number.denominator_.toString();
  }
  return stringRepresentation;
}

string Rational::asDecimal(size_t precision) const {
  string stringRepresentation;
  if (!isNonNegative_) {
    stringRepresentation += "-";
  }
  BigInteger base = numerator_;
  base.restoreDegree(precision / kBaseModule);
  base *= std::pow(10, precision % kBaseModule);

  string temp = (base / denominator_).toString();
  if ((temp.length()) < precision) {
    stringRepresentation += "0.";
    char* zeros = new char[precision - temp.length()];
    std::memset(zeros, '0', precision - temp.length());
    zeros[precision - temp.length()] = '\0';
    stringRepresentation += zeros;
    stringRepresentation += temp;
    delete[] zeros;
  } else {
    stringRepresentation += temp.substr(0, temp.length() - precision);
    if (bool(denominator_) == 1 && precision != 0) {
      stringRepresentation += ".";
      stringRepresentation += temp.substr(temp.length() - precision, precision);
    }
  }
  return stringRepresentation;
}

std::ostream& operator<<(std::ostream& OutputStream,
                         const BigInteger& integer) {
  OutputStream << integer.toString();
  return OutputStream;
}

std::istream& operator>>(std::istream& InputStream, BigInteger& integer) {
  std::string str;
  InputStream >> str;
  integer = BigInteger(str);
  return InputStream;
}

BigInteger operator""_bi(const char* str, size_t) { return BigInteger(str); }
