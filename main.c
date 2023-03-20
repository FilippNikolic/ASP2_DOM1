#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <time.h>
#ifndef STACK_H
#define STACK_H

typedef struct st
{
    void **content;
    int top;
    int max_size;
} Stack;

Stack *createStack(int initial_size);
void emptyStack(Stack *);
void eraseStack(Stack *);
int push(Stack *, void *);
void *pop(Stack *);
int isEmptyStack(Stack *);
#endif

int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

Stack *createStack(int initial_size)
{
    if( initial_size < 1 ) return 0;
    else
    {
        Stack *ret = (Stack *)calloc(1, sizeof(Stack));

        if( ! ret ) return 0;
        ret->content = calloc(initial_size, sizeof(void *));
        if( ! ret->content )
        {
            free(ret);
            return 0;
        }
        ret->max_size = initial_size;
        return ret;
    }
}

void emptyStack(Stack *s)
{
    if( s )  s->top = 0;
}

void eraseStack(Stack *s)
{
    if( s )
    {
        free(s->content);
        free(s);
    }
}

int push(Stack *s, void *dat)
{
    if( ! s ) return 0;

    if( s->top == s->max_size*9/10 )
    {
        Stack **tmp=(Stack **)realloc(s->content,sizeof(void *)*s->max_size*2);
        if( ! tmp ) return 0;

        s->content = tmp;
        s->max_size *= 2;
    }
    s->content[s->top++] = dat;
    return 1;
}

void *pop(Stack *s)
{
    if( ! s ) return 0;
    if( s->top == 0 ) return 0;
    return s->content[--s->top];
}

int isEmptyStack(Stack *s)
{
    if( ! s )  return 1;
    return s->top == 0;
}

#ifndef QUEUE_H
#define QUEUE_H
typedef struct queue
{
    void **content;
    int front, rear;
    int max_size;
} Queue;

Queue *createQueue(int initial_size);
void emptyQueue(Queue *);
void eraseQueue(Queue *);
int insertIntoQueue(Queue *, void *);
void *removeFromQueue(Queue *);
int isEmptyQueue(Queue *);
#endif

Queue *createQueue(int initial_size)
{
    Queue *ret = (Queue *)calloc(1, sizeof(Queue));

    if( ! ret ) return 0;
    ret->content = calloc(initial_size, sizeof(void *) );
    if( ! ret->content )
    {
        free(ret);
        return 0;
    }
    ret->max_size = initial_size;
    return ret;
}

void emptyQueue(Queue *q)
{
    if( q ) q->front = q->rear = 0;
}

void eraseQueue(Queue *q)
{
    if( q )
    {
        free(q->content);
        free(q);
    }
}

int insertIntoQueue(Queue *q, void *dat)
{
    if( ! q ) return 0;
    if( (q->rear + 1)%q->max_size == q->front )
    {
        void **c = (void **)calloc(q->max_size*2, sizeof(void *));
        int i, j;

        if( ! c )   return 0;
        for(i = 0, j = q->front; j != q->rear; j = (j+1) % q->max_size, i++)
            c[i] = q->content[j];

        q->front = 0; q->rear = i;
        free(q->content);
        q->content = c;
        q->max_size *= 2;
    }
    q->content[q->rear++] = dat;
    q->rear = q->rear % q->max_size;
    return 1;
}

void *removeFromQueue(Queue *q)
{
    void *ret = 0;
    if( ! q ) return 0;
    if( q->front == q->rear ) return 0;
    ret = q->content[q->front++];
    q->front = q->front % q->max_size;
    return ret;
}

int isEmptyQueue(Queue *q)
{
    if( ! q )  return 1;
    return q->front == q->rear;
}

#ifndef BINTREE_H
#define BINTREE_H
typedef struct btNode
{
    int key,vis;
    struct btNode *parent, *left, *right;
} BinTreeNode;

typedef struct binTree
{
    BinTreeNode *root;
    int treeHeight;
} BinTree;

typedef struct btVisitor
{
    BinTree     *tree;
    Stack       *stack;
    BinTreeNode *current;
} BinTreeVisitor;

/* tree creation, emptying and deletion */
BinTree *createBinTree();
void deleteBinTree(BinTree *);
void emptyBinTree(BinTree *);

/* key (node) lookup:
- node equals NULL if key is not found
- parent is the parent node of the node where the key is or should be
returns the level at which the key has been found  */
int findKey(BinTree *, int, BinTreeNode **, BinTreeNode **);

/* key insertion and deletion */
int insertKey(BinTree *, int);
int deleteKey(BinTree *, int);

/* tree visiting */
BinTreeVisitor *createVisitor(BinTree *);
void deleteVisitor(BinTreeVisitor *);
void visitPreOrder(BinTreeVisitor *);
void visitInOrder(BinTreeVisitor *);
void visitPostOrder(BinTreeVisitor *);
BinTreeNode *successor(BinTreeNode *node);

/* tree display */
void printTree(BinTree *tree);
void printVisit(BinTree *, void (*visitMethod)(BinTreeVisitor *) );

#endif

BinTree *createBinTree() {
    BinTree *bt = (BinTree *)calloc(1, sizeof(BinTree) );
    if( ! bt )  return 0;
    bt->treeHeight = -1;
    return bt;
}

void deleteBinTree(BinTree *bt) {
    emptyBinTree(bt);
    free(bt);
}

void emptyBinTree(BinTree *bt){
    if( bt ) {
        BinTreeVisitor *visitor = createVisitor(bt);
        if( ! visitor )   return;
        visitPostOrder(visitor);
        while( visitor->current ) {
            free(visitor->current);
            visitPostOrder(visitor);
        }
        bt->root = 0;
        bt->treeHeight = -1;
    }
}

int subtreeHeight(BinTreeNode *node) {
    if( ! node )      return -1;
    return 1+max( subtreeHeight(node->left),subtreeHeight(node->right));
}

int treeHeight(BinTree *bt) {
    if( ! bt ) return -1;
    if( ! bt->root ) return -1;
    return subtreeHeight(bt->root);
}

int findKey(BinTree *tree,int key,BinTreeNode **node,BinTreeNode **parent)
{
    int level = 0;
    *node = *parent = 0;
    if( ! tree )      return -1;

    *node = tree->root;
    while(*node)
    {
        if((*node)->key == key) break;
        level++;
        *parent = *node;
        if( (*node)->key > key ) *node = (*node)->left;
        else                     *node = (*node)->right;
    }
    return level;
}

int insertKey(BinTree *tree, int key)
{
    if( ! tree )  return 0;
    else {
        BinTreeNode *node, *parent, *newnode;
        int level = findKey(tree, key, &node, &parent);

        if( node ) {
            parent = successor(node);
            if(parent == NULL){
                parent = node;
            }
        }
        newnode = calloc(1, sizeof(BinTreeNode) );
        if( !newnode ) return 0;

        newnode->key = key;
        newnode->parent = parent;

        if( parent ){
            if( parent->key > key ) parent->left = newnode;
            else                    parent->right = newnode;
        }
        else                        tree->root = newnode;

        tree->treeHeight = max(tree->treeHeight, level);

        return 1;
    }
}

BinTreeNode *predecessor(BinTreeNode *node)
{
    if( ! node ) return 0;
    node = node->left;
    while( node && node->right ) node=node->right;
    return node;
}

BinTreeNode *successor(BinTreeNode *node)
{
    if( ! node ) return 0;
    node = node->right;
    while( node && node->left ) node = node->left;
    return node;
}

void linkParent(BinTreeNode *parent, BinTreeNode *node, BinTreeNode *next)
{
    if( parent && node )
    {
        if( parent->left == node )  parent->left = next;
        else                        parent->right = next;
    }
}

int deleteKey(BinTree *tree, int key)
{
    if( ! tree )    return 0;
    else
    {
        BinTreeNode *node, *parent, *succ;
        int level = findKey(tree, key, &node, &parent);

        if( ! node )      return 0;

        if( ! node->left && ! node->right )
        {
            if( parent )  linkParent(parent, node, 0);
            else          tree->root = 0;
            free(node);
        }
        else  if( ! (node->left && node->right) )
        {
            BinTreeNode *next;

            if( node->left )  next = node->left;
            else              next = node->right;

            if( parent )      linkParent( parent, node, next );
            else              tree->root = next;

            return 1;
        }
        else  {
            int tmp;
            succ = successor(node);
            tmp = succ->key;
            succ->key = node->key;
            node->key = tmp;

            if( succ->parent != node ) succ->parent->left = succ->right;
            else                       node->right = succ->right;
            free(succ);
        }
        tree->treeHeight = treeHeight(tree);
        return 1;
    }
}

/* tree display */
void printTree(BinTree *tree) {
    if( ! tree )   return;
    if( ! tree->root )  return;
    else {
        Queue *q = createQueue( tree->treeHeight*tree->treeHeight);
        int i, line_len = 62;
        int first_skip = line_len, in_between_skip;

        if( ! q )   return;
        insertIntoQueue(q, tree->root);
        for( i = 0; i <= tree->treeHeight; i++ ) {
            int j = 1 << i, k, l;
            in_between_skip = first_skip;
            first_skip = (first_skip-2)/2;
            for( k = 0; k < first_skip; k++) putchar(' ');
            for(k = 0; k < j; k++) {
                BinTreeNode *btn = (BinTreeNode *)removeFromQueue(q);
                if( btn ) {
                    insertIntoQueue(q, btn->left);
                    insertIntoQueue(q, btn->right);
                } else  {
                    insertIntoQueue(q, 0);
                    insertIntoQueue(q, 0);
                }
                if( btn )  printf("%2d", btn->key );
                else       printf("  ");
                for( l = 0; l < in_between_skip; l++) putchar(' ');
            }
            putchar('\n');
            putchar('\n');
        }
        eraseQueue(q);
    }
}

void printVisit(BinTree *bt, void (*visitMethod)(BinTreeVisitor *) ) {
    BinTreeVisitor *v = createVisitor(bt);
    if( ! v )  return;
    do {
        visitMethod(v);
        if( v->current ) printf("%3d", v->current->key );
    } while( v->current );
    deleteVisitor(v);
    printf("\n");
}

BinTreeVisitor *createVisitor(BinTree *bt) {
    if( ! bt ) return 0;
    else {
        BinTreeVisitor *ret=(BinTreeVisitor*)calloc(1,sizeof(BinTreeVisitor));
        if( ! ret ) return 0;

        ret->tree = bt;
        ret->stack = createStack(bt->treeHeight*2);
        if( ! ret->stack ) {
            free(ret);
            return 0;
        }
        return ret;
    }
}

void deleteVisitor(BinTreeVisitor *v) {
    if( v ) {
        free(v->stack);
        free(v);
    }
}

void visitPreOrder(BinTreeVisitor *v) {
    if( v )  {
        BinTreeNode *next;
        if( isEmptyStack(v->stack) && v->current && !v->current->left) {
            v->current = 0;
            return;
        }
        if( ! v->current ) push(v->stack, v->tree->root);
        next = (BinTreeNode *)pop(v->stack);
        v->current = next;
        if( next ) {
            if( next->right ) push(v->stack, next->right);
            if( next->left )  push(v->stack, next->left);
        }
    }
}

void visitInOrder(BinTreeVisitor *v) {
    if( v ) {
        BinTreeNode *next;

        if( isEmptyStack(v->stack) && v->current ) {
            v->current = 0;
            return;
        }

        if( ! v->current ) {
            next = v->tree->root;
            while( next ) {
                push(v->stack, next);
                next = next->left;
            }
        }

        next = (BinTreeNode *)pop(v->stack);
        v->current = next;

        next = next->right;
        while(next){
            push(v->stack, next);
            next = next->left;
        }
    }
}

void visitPostOrder(BinTreeVisitor *v) {
    if( v ) {
        BinTreeNode *next;

        if( isEmptyStack(v->stack) && v->current ) {
            v->current = 0;
            return;
        }

        if( ! v->current ) {
            next = v->tree->root;
            while( next ) {
                push(v->stack, next);
                push(v->stack, next);
                next = next->left;
            }
        }

        while( ! isEmptyStack(v->stack) )
        {
            next = (BinTreeNode *)pop(v->stack);
            if( next ) {
                push(v->stack, 0);
                next = next->right;
                while(next)  {
                    push(v->stack, next);
                    push(v->stack, next);
                    next = next->left;
                }
            }  else {
                next = (BinTreeNode *)pop(v->stack);
                v->current = next;
                return;
            }
        }
        v->current = 0;
    }
}

int m_pretraga(int m,int key,int niz[],int k){
    int l = 0,mid1,mid2,high = k - 1,low = 0;
    mid1 = (low+high)/m;
    if(key < niz[0] || key > niz[k-1] )return -1;
    while (low <= high){
        if(high == low){
            if(niz[high]==key) {
                return high + 1;
            } else{
                return -1;
            }
        }
        mid2 = mid1;
        for (int i = 0; i < m; ++i) {
            if(key == niz[mid2 + l]){
                return mid2+l+1;
            }
            if(key < niz[mid2+l]){
                high = mid2+l-1;
                mid1 = (high - low)/m;
                if(mid1 == 0){
                    if(niz[low] == key){
                        return low+1;
                    }
                    mid1 = 1;
                }
                l = low;
                break;
            } else{
                low = (mid2+l) + 1;
                if(i == m - 2){
                    mid2 = high-l;
                    continue;
                }
                mid2  += mid1;
            }
        }
    }
    return -1;
}

int m_pretraga_opt(int m,int key,int niz[],int k){
    int l = 0,mid1,mid2,high = k - 1,low = 0;
    mid1 = low+(high-low)*(key - niz[low])/(niz[high]-niz[low]);
    if(key < niz[0] || key > niz[k-1] )return -1;
    while (low <= high){
        if(high == low){
            if(niz[high]==key) {
                return high + 1;
            } else{
                return -1;
            }
        }
        mid2 = mid1;
        for (int i = 0; i < m; ++i) {
            if(key == niz[mid2 + l]){
                return mid2+l+1;
            }
            if(key < niz[mid2+l]){
                high = mid2+l-1;
                mid1 = (high - low)/m;
                if(mid1 == 0){
                    if(niz[low] == key){
                        return low+1;
                    }
                    mid1 = 1;
                }
                l = low;
                break;
            } else{
                low = (mid2+l) + 1;
                if(i == m - 2){
                    mid2 = high-l;
                    continue;
                }
                mid2  += mid1;
            }
        }
    }
    return -1;
}

void obrada_rez(int rez,int niz[],int tmp[]){
    while ((niz[rez] == niz[rez-1]) && rez > 1){
        rez--;
    }
    while ((niz[rez-1] == niz[rez]) && (tmp[rez-1]==1)){
        rez++;
    }
    if(niz[rez-1] != niz[rez] && tmp[rez-1]==1){
        rez = -1;
    }
    if(rez == -1){
        printf("Neuspesna pretraga za broj %d!\n",niz[rez-1]);
    } else{
        printf("%d -> %d ",niz[rez-1],rez);
        tmp[rez-1] = 1;
        printf("\n");
    }

}

float loga_b(float a, float b){
    float t;
    t = log2f(a) / log2f(b);
    return t;
}


int main() {

    int zad = 1, sdr, k, obs,low = 0,high,mid;
    printf("Duzina niza je:");
    scanf("%d", &k);
    int niz[k];
    switch (zad) {
        case 1: {
            //zadatak 1
            while (1){
                if(sdr == 7){
                    break;
                }
                printf("1. Unos uredjenog niza\n"
                       "2. Generisanje uredjenog niza\n"
                       "3. Pretraga uredjenog niza\n"
                       "4. Pretraga niza kljuceva\n"
                       "5. Evaluacija preformanse m-arnog pretrazivanje\n"
                       "6. Evaluacija preformanse optimizovanog m-arnog pretrazivanje\n"
                       "7. Prelazak na zadatak 2.\n"
                       "8. Kraj programa\n");
                scanf("%d", &sdr);
                switch (sdr) {
                    case 1: {
                        printf("Unos niza iz datoteke(1) ili standardnog ulaza(2):");
                        int p;
                        scanf("%d", &p);
                        if (p != 2 && p != 1) {
                            perror(NULL);
                        }
                        if (p == 1) {
                            FILE *file;
                            file = fopen("C:\\Users\\nfili\\OneDrive\\Desktop\\ASP\\ASP2\\K1\\Dom_Asp2_1\\dat.txt","r");
                            int t,dat[100],br = 0;
                            do{
                                fscanf(file,"%d, ",&t);
                                dat[br++] = t;
                            } while (br<100);
                            fclose(file);
                            for (int i = 0; i < k-1; ++i) {
                                niz[i] = dat[i];
                                printf("%d->",niz[i]);
                            }
                            printf("%d\n",dat[k-1]);
                        }
                        if (p == 2) {
                            printf("Unesite niz:");
                            for (int i = 0; i < k; ++i) {
                                scanf("%d", &niz[i]);
                            }
                        }
                        break;
                    }
                    case 2: {
                        int top,bottom;
                        printf("Unesite obsek niza(od - do):");
                        scanf("%d %d", &bottom,&top);
                        srand(time(0));
                        for (int i = 0; i < k; ++i) {
                            niz[i] = (rand() % (top - bottom + 1)) + bottom;
                        }
                        for (int i = 0; i < k; ++i) {
                            for (int j = i + 1; j < k; ++j) {
                                if (niz[i] > niz[j]) {
                                    int tmp = niz[i];
                                    niz[i] = niz[j];
                                    niz[j] = tmp;

                                }
                            }
                        }
                        for (int i = 0; i < k-1; ++i) {
                            printf("%d->", niz[i]);
                        }
                        printf("%d\n\n",niz[k-1]);
                        break;
                    }
                    case 3: {
                        int l = 0,flag = 0,tmp[k];
                        printf("Unesite m i kljuc:");
                        int m,key,rez;
                        scanf("%d %d",&m,&key);
                        for (int i = 0; i < k; ++i) {
                            tmp[i] = 0;
                        }
                        high = k - 1;
                        mid = (low+high)/m;
                        rez = m_pretraga(m,key,niz,k);
                        obrada_rez(rez,niz,tmp);
                        break;
                    }
                    case 4:{
                        int p,m,rez;
                        printf("Kolika je velicina niza kljuceva:");
                        scanf("%d",&p);
                        printf("Kolika je m:");
                        scanf("%d",&m);
                        int P[p],I[p];
                        printf("Unesite niz kljuceva:");
                        int tmp[k];
                        for (int i = 0; i < k; ++i) {
                            tmp[i] = 0;
                        }
                        for (int i = 0; i < p; ++i) {
                            scanf("%d",&P[i]);
                            rez = m_pretraga(m,P[i],niz,k);
                            I[i] = rez;
                        }
                        for (int i = 0; i < p; ++i) {
                            obrada_rez(I[i],niz,tmp);
                        }
                        printf("\n");
                        break;
                    }
                    case 5: {
                        int p,rez,tmp[k];
                        for (int i = 0; i < k; ++i) {
                            tmp[i] = 0;
                        }
                        printf("Uneste duzinu niza kljuceva:");
                        scanf("%d", &p);
                        int tek[p];
                        for (int i = 0; i < p; ++i) {
                            printf("Unesite vrednost kljuca:");
                            scanf("%d", &tek[i]);
                            rez = m_pretraga(2,tek[i],niz,k);
                            obrada_rez(rez,niz,tmp);
                            float r = k;
                            for (int j = 2; j < 7; ++j) {
                                float l = j;
                                float ev = loga_b(r,l);
                                printf("Za parametar vrednosti %d prosecan broj koraka je %f\n",j,ev);
                            }
                        }
                        printf("\n");
                        break;
                    }
                    case 6: {
                        int p, rez, tmp[k],m;
                        for (int i = 0; i < k; ++i) {
                            tmp[i] = 0;
                        }
                        printf("Uneste duzinu niza kljuceva:");
                        scanf("%d", &p);
                        int tek[p];
                        int I[p];
                        for (int i = 0; i < k; ++i) {
                            tmp[i] = 0;
                        }
                        printf("Unesite vrednosti kljuceva:");
                        for (int i = 0; i < p; ++i) {
                            scanf("%d", &tek[i]);
                            rez = m_pretraga(2,tek[i],niz,k);
                            I[i] = rez;
                        }
                        float l,r;
                        int tmp2 = k;
                        for (int i = 0; i < p; ++i) {
                            obrada_rez(I[i],niz,tmp);
                        }
                        for (int z = 2; z <= 6; ++z){
                            float ev = 0;
                            for (int i = 0; i < p; ++i) {
                                r = tmp2;
                                l = z;
                                ev += loga_b(r, l);
                                if(I[i] != -1){
                                    tmp2 = k - I[i];
                                }
                            }
                            tmp2 = k;
                            ev = ev/p;
                            printf("Za parametar vrednosti %d prosecan broj koraka je %f\n",z, ev);
                        }
                        printf("\n");
                        break;
                    }
                    case 7:{
                        break;
                    }
                    default:{
                        return 0;
                    }
                }

            }
        }
        case 2: {
            //zadatak2
            int nmb,niz2[100],cnt=0;
            for (int i = 0; i < 100; ++i) {
                niz2[i] = 0;
            }
            BinTree *tree;
            while (1){
                printf("1. Formiranje BST-a\n"
                       "2. Umetanje novog kljuca u stablo\n"
                       "3. Pretraga stabla na zadati kljuc\n"
                       "4. Ispis stabla\n"
                       "5. Prosecan broj koraka za pretrazivanje u zadatom stablu\n"
                       "6. Brisanje stabla iz memorije\n"
                       "7. Kraj programa\n");
                scanf("%d", &nmb);
                switch (nmb) {
                    case 1: {
                        tree = createBinTree();
                        for (int i = 0; i < k; ++i) {
                            insertKey(tree, niz[i]);
                        }
                        break;
                    }
                    case 2: {
                        int tmp;
                        printf("Unesite kljuc koji zelite da umetnete u stablo:");
                        scanf("%d",&tmp);
                        insertKey(tree,tmp);
                        niz2[cnt++] = tmp;
                        break;
                    }
                    case 3:{
                        int tmp;
                        printf("Unesite kljuc koji zelite da nadjete u stablu:");
                        scanf("%d",&tmp);
                        struct btNode* tek = tree->root;
                        while (tek){
                            if(tek->key==tmp){
                                printf("%d\n",tek->key);
                                break;
                            }
                            printf("%d->",tek->key);
                            if(tek->key>tmp){
                                tek = tek->left;
                            } else{
                                tek = tek->right;
                            }
                        }
                        if(!tek){
                            printf("NULL\n");
                        }
                        break;
                    }
                    case 4:{
                        printTree(tree);
                        break;
                    }
                    case 5:{
                        float br=0;
                        struct btNode* tek = tree->root;
                        for (int i = 0; i < k; ++i){
                            while (tek) {
                                br++;
                                if (tek->key == niz[i]) {
                                    /*while(tek->vis == 1){
                                        tek = tek->right;
                                        br++;
                                        while( tek->key != niz[i] ) { tek = tek->left;br++; }
                                    }
                                    tek->vis = 1;*/
                                    break;
                                }
                                br++;
                                if (tek->key > niz[i]) {
                                    tek = tek->left;
                                } else {
                                    tek = tek->right;
                                }
                            }
                            tek = tree->root;
                        }
                        tek = tree->root;
                        for (int i = 0; i < cnt; ++i){
                            while (tek) {
                                br++;
                                if (tek->key == niz2[i]) {
                                    /*while(tek->vis == 1){
                                        tek = tek->right;
                                        br++;
                                        while( tek->key != niz2[i] ) { tek = tek->left;br++; }
                                    }
                                    tek->vis = 1;*/
                                    break;
                                }
                                br++;

                                if (tek->key > niz2[i]) {
                                    tek = tek->left;
                                } else {
                                    tek = tek->right;
                                }
                            }
                            tek = tree->root;
                        }
                        float opt;
                        opt = br/(k+cnt);
                        printf("Prosecan broj koraka za nalazenje kljuca u zadatom stablu je %f\n\n",opt);
                        break;
                    }
                    case 6:{
                        deleteBinTree(tree);
                        return 0;
                    }
                    default:
                       return 0;
                }
            }

        }
    }
}