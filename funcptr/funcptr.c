#include <stdio.h>
#include <stdlib.h>

typedef struct _Animal Animal;

struct _Animal{
    int age;
    char *name;
    void (*funcptr)(Animal *);   
};

void Show_info(Animal *animal)
{
    printf("%s < I am %d years old!!\n", animal->name, animal->age);
}

void walking(Animal *animal)
{
    printf("%s < I am walking!!\n", animal->name);
}

int main(void) 
{   
    Animal *animal = malloc(3 * sizeof (Animal));
    animal[0].age = 7;
    animal[0].name = "mona";
    animal[0].funcptr = Show_info;

    animal[1].age = 11;
    animal[1].name = "chiana";
    animal[1].funcptr = walking;

    animal[2].funcptr = NULL;

    printf("Pattern 1\n");
    for(Animal *p = animal; p->funcptr; p+=1)
    {
        p->funcptr(p);
    }

    printf("Pattern 2\n");
    for(int i=0;i<2;i++)
    {
        animal[i].funcptr(&animal[i]);
    }
    return 0;
}