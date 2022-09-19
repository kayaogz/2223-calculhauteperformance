#include <stdio.h>
#include <unistd.h>
#include "omp.h"

int comptwo()
{
  sleep(2);
  return 2;
}

int compthree()
{
  sleep(3);
  return 3;
}

int main()
{
  int two, three, five = 0;

#pragma omp parallel num_threads(4)
  {
    int thid = omp_get_thread_num();
    int numth = omp_get_num_threads();
    printf("hello from %d/%d\n", thid, numth);
    //
#pragma omp sections
    {
#pragma omp section
      {
        printf("section two th %d\n", thid);
        two = comptwo();
        five += two;

      }
#pragma omp section
      {
        printf("section three th %d\n", thid);
        three = compthree();
//#pragma omp atomic
#pragma omp critical
        {
          five += three;
        }
      }
#pragma omp section
      {
        printf("asdfasdf\n");
      }
    }
  }
  five = two + three;
  printf("five = %d\n", five);

  return 0;
}
