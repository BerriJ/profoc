//  Student.h

#ifndef conline_h
#define conline_h

#include <string>
#include <vector>

class Student
{
public:
  // Constructor
  Student(std::string name, int age, bool male);

  // Getters
  std::string GetName();
  int GetAge();
  bool IsMale();
  std::vector<int> GetFavoriteNumbers();

  // Methods
  bool LikesBlue();

private:
  // Member variables
  std::string name;
  int age;
  bool male;
  std::vector<int> favoriteNumbers;
};

#endif /* conline_h */
