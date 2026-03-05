#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* This program simulates the game of craps. */
/* Intructions:*/
/* 1. Roll two dice and calculate their sum. */
/* 2. If the sum is 7 or 11 on the first throw, you win. */
/* 3. If the sum is 2, 3, or 12 on the first throw, you lose. */
/* 4. If the sum is 4, 5, 6, 8, 9, or 10 on the first throw, 
then that sum becomes your "point". */
/* 5. To win, you must continue rolling the dice until 
you "make your point" (i.e., roll that same point value). */ 


/* Define clases and functions */

enum Status{ WIN, LOSE, CONTINUE};
int rollDice();
int calculateStatus(int sum, int point);

int main() {

    /* Modify random seed */
    srand((unsigned)time(NULL));

    /* Define variables */
    enum Status game_status ;
    int obtained_number;
    int point = 0;

    /* Roll dice*/
    printf("Presiona enter para lanzar los dados...");
    getchar();
    obtained_number = rollDice();

    game_status = calculateStatus(obtained_number, point);


    point = obtained_number;
    printf("Your point is %d\n", point);

    while (game_status == CONTINUE) {
        
        printf("Presiona enter para lanzar los dados...");
        getchar();
        obtained_number = rollDice();
        printf("Your point is %d\n", point);

        if (obtained_number == point) {
            game_status = WIN;
            printf("You win!\n");
        } else if (obtained_number == 7) {
            game_status = LOSE;
            printf("You lose!\n");
        }
    }

    return 0;
}


int rollDice() {

    int die1 = 1 + rand() % 6;
    int die2 = 1 + rand() % 6;

    int sum = die1 + die2;
    printf("You rolled %d + %d = %d\n", die1, die2, sum);
    return sum;
}

int calculateStatus(int sum, int point) {

    if (sum == 7 || sum == 11) {
        printf("You win!\n");
        exit(0);
        return WIN;
    } else if (sum == 2 || sum == 3 || sum == 12) {
        printf("You lose!\n");
        exit(0);
        return LOSE;
    } else {
        return CONTINUE;
    }
}