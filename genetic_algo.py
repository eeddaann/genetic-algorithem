import numpy as np
import matplotlib.pyplot as plt

width=10
n=10
best_lst=[]
best_sol=0
def f(x):
    '''
    calc objerctive function
    :param x: a tuple of x1.x2 coordinates
    :return: the value of objective function
    '''
    if x[0]==0 or x[1]==0:
        return 0.000001 #prevents log(0) instead returns a small number that probably won't be the optimal solution

    if x[0]+x[1]>3000 or x[0]<300 or (float(x[0])/x[1])<1.5: #constraints
        return 0.000001
    return 0.7*np.log(x[0]/300)+0.3*np.log(x[1]) #calc of objective function

def f_chrom(chrom):
    '''
    :param chrom: a chromosome
    :return: objective function value of the chromosome
    '''
    gene1,gene2=chrom_to_gene(chrom)
    x=gene_to_x(gene1),gene_to_x(gene2)
    return f(x)


def x_to_gene(x):
    '''
    :param x: an integer
    :return: the representation of x as a gene
    '''
    return str(np.binary_repr(x,width=width))

def gene_to_x(gene):
    '''
    :param gene: a gene
    :return: the original integer
    '''
    gene=str(gene)
    return int(gene,2)

def gene_to_chrom(g1,g2):
    '''
    :param g1: a gene
    :param g2: a gene
    :return: a chromosome
    '''
    return g1+g2

def chrom_to_gene(chrom):
    '''
    :param chrom: a chromosome
    :return: a tuple of the genes
    '''
    gene1=chrom[0:10]
    gene2=chrom[10:]
    return gene1,gene2

def random_chrom():
    '''
    generates a random chromosome
    '''
    chrom=''
    for i in range(2*width):
        chrom+= np.random.choice(['0','1'])
    return chrom

def init_pop(n=10):
    '''
    :param n: the population size
    :return: a list of chromosomes
    '''
    return [random_chrom() for i in range(n)]

def calc_fit(pop):
    val_func=[]
    for i in pop:
        val_func.append(f_chrom(i))# calc fit of each chromosome in population
    sum_func=sum(val_func)#clac sum of values
    fit_vals=[]
    for j in val_func:
        fit_vals.append(float(j)/sum_func)# divide each chromosome to get proportion
    return fit_vals


def mutate(chrom,p=0.05):
    '''
    :param chrom: a chromosome
    :param p: the probability for a mutation in each bit
    :return:a mutated chromosome
    '''
    chrom=list(chrom)
    for i in range(len(chrom)):
        j=np.random.rand(0,1)
        if j<=p:
            chrom[i]=1-chrom[i]
    return "".join(chrom)

def twin_sperator(chrom):
    '''
    to assure a diverse population this function mutate a random bit in the given chromosome
    :param chrom: a chromosome
    :return: a mutated chromosome
    '''
    chrom=list(chrom)
    point=np.random.choice(range(len(chrom)-1))
    chrom[point]=str(1-int(chrom[point]))
    return "".join(chrom)

def select_parents(pop):
    '''
    :param pop: list of chromosomes
    :return: a list of tuples that represents parents
    '''
    parents=[]
    a=calc_fit(pop)#the proportions
    for i in range(len(pop)/2):
        par1=np.random.choice(pop,p=a)
        par2=par1 #to prevent same chromosome crossover
        while par2==par1:
            par2=np.random.choice(pop)
        parents.append((par1,par2))
    return parents

def crossover(parents):
    point=np.random.randint(1,(2*n)-1) #get a crossover point
    child1=parents[0][0:point]+parents[1][point:]
    child2=parents[0][point:]+parents[1][0:point]
    return [child1,child2]

def offsprings(parents):
    children=[]
    for i in parents:
        cur_children=crossover(i)
        while cur_children[0]  in children or cur_children[1] in children:# check if one or more of the children are already in the children list
            if cur_children[0]  in children:
                cur_children[0]=twin_sperator(cur_children[0])
            if cur_children[1]  in children:
                cur_children[1]=twin_sperator(cur_children[1])
        children.extend(cur_children)
    return children

def survival(pop,children):
    pop_val=[f_chrom(p) for p in pop] #objective function value of each chromosome in population
    children_val=[f_chrom(ch) for ch in children] #objective function value of each chromosome in children
    best_parent_val=max(pop_val) # in population of ten a ten precent is exactly one, so only one parents will move to the next generation (elitism)
    best_parent =pop[pop_val.index(best_parent_val)]
    worst_child_val = min(children_val)
    children[children_val.index(worst_child_val)]=best_parent
    children_val[children_val.index(worst_child_val)]=best_parent_val
    best_lst.append(max(children_val))
    global best_sol
    if best_lst[-1]!=best_sol:
        best_sol=children[children_val.index(best_lst[-1])]
    return children




pop=init_pop()
for gen in range(1000):#1000 generations
    parents=select_parents(pop)
    offspring=offsprings(parents)
    mutated_offspring=[mutate(off) for off in offspring]#mutate offsprings
    pop=survival(pop, mutated_offspring)#the new population is the chromosomes that survived


#create the graph
x=range(1000)
y=best_lst
plt.plot(x, y)
plt.show()
print "best solution found at:",[gene_to_x(i)for i in chrom_to_gene(best_sol)]#print best solution
print "with objective function value of:",best_lst[-1]


