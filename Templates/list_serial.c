#include <stdio.h>
#include <stdlib.h>

#define LIST_ELEMS 20

int count = 0;

typedef struct node node;

struct node {
  int data;
  int id;
  node *next;
};

void initialize(node *p){
  p->data = count;
  count++;
}

void process(node *p){
  p->id = 1;
}

int main(void){

  node *head;
  node *p;
  int i;

  /* create list */
  head = (node *) malloc (sizeof(node));
  p = head;
  for(i=0; i<LIST_ELEMS; ++i){
    p->next = (node *) malloc (sizeof(node));
    p = p->next;
  }
  p->next = NULL;

  p = head;
  /* initialize list items */
  while (p) {
    initialize(p);
    p = p->next;
  }  

  p = head;
  /* process list items */
  while (p) {
    process(p);
    p = p->next;
  }  

  p = head;
  /* print list */
  while (p) {
    printf("num: %d  -> %d  \n",p->data, p->id);
    p = p->next;
  }
  printf("\n");

  return 0;
}
