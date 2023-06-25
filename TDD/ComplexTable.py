"""Define global variables"""

epi=1e-6

complex_table = dict()

complex_entry_table =dict()

cacheCount = 1


class ComplexTableEntry:
    def __init__(self,val=0):
        self.val = val
        
    def __str__(self):
        return str(val)
        
class Complex:
    def __init__(self,c=0):
        self.r = ComplexTableEntry(c.real)
        self.i = ComplexTableEntry(c.imag)
        self.next = None
        
    def __add__(self,other):
        res = cacheAvail
        res.r.val = self.r.val+other.r.val
        res.i.val = self.i.val+other.i.val
        return res
    
    def __sub__(self,other):
        res = cacheAvail
        res.r.val = self.r.val-other.r.val
        res.i.val = self.i.val-other.i.val
        return res    
    
    def __mul__(self,other):
        
        if self==cn0 or other == cn0:
            return cn0         
        
        res = cacheAvail
        if self==cn1:
            res.r.val = other.r.val
            res.i.val = other.i.val
            return res
        
        if other==cn1:
            res.r.val = self.r.val
            res.i.val = self.i.val
            return res
    
    
        ar=self.r.val
        ai=self.i.val
        br=other.r.val
        bi=other.i.val
        
        res.r.val = ar*br-ai*bi
        res.i.val = ar*bi+ai*br
        return res
    
    def __truediv__(self,other):
        
        if self==other:
            return cn1
        
        if self == cn0:
            return cn0
        
        res = cacheAvail
        
        if other==cn1:
            res.r.val = self.r.val
            res.i.val = self.i.val
            return res        
        
        ar=self.r.val
        ai=self.i.val
        br=other.r.val
        bi=other.i.val
        
        cmag = br*br+bi*bi
        res.r.val = (ar*br+ai*bi)/cmag
        res.i.val = (ai*br-ar*bi)/cmag
        return res   
    
    def norm(self):
        ar=self.r.val
        ai=self.i.val
        return ar*ar+ai*ai
    
    
    def __eq__(self,other):
        if self.r == other.r and self.i == other.i:
#             print(True,self,other,id(self.r),id(other.r),id(self)==id(other))
            return True
        else:
#             print(False,self,other,id(self.r),id(other.r))
            return False
        
    def __str__(self):
        return str(self.r.val+1j*self.i.val)        
        
cn0 = Complex(0)
cn1 = Complex(1)
cacheAvail = Complex()

def cn_mul(res:Complex, a:Complex, b:Complex):
    """res=a*b"""
    if equalsOne(a):
        res.r.val = b.r.val
        res.i.val = b.i.val
        return
    if equalsOne(b):
        res.r.val = a.r.val
        res.i.val = a.i.val
        return 
    if equalsZero(a) or equalsZero(b):
        res.r.val = 0
        res.i.val = 0
        return 

    ar=a.r.val
    ai=a.i.val
    br=b.r.val
    bi=b.i.val
        
    res.r.val = ar*br-ai*bi
    res.i.val = ar*bi+ai*br

def cn_mulCached(a:Complex,b:Complex):
#     res = getCachedComplex()
    global cacheAvail,cacheCount
    res = cacheAvail
    cacheAvail=cacheAvail.next
    cacheCount-=1
    if equalsOne(a):
        res.r.val = b.r.val
        res.i.val = b.i.val
        return res
    if equalsOne(b):
        res.r.val = a.r.val
        res.i.val = a.i.val
        return res
    if equalsZero(a) or equalsZero(b):
        res.r.val = 0
        res.i.val = 0
        return res
    ar=a.r.val
    ai=a.i.val
    br=b.r.val
    bi=b.i.val
    res.r.val = ar*br-ai*bi
    res.i.val = ar*bi+ai*br
    return res

def cn_add(res:Complex, a:Complex, b:Complex):
    """res=a*b"""
    res.r.val = a.r.val+b.r.val
    res.i.val = a.i.val+b.i.val

def cn_addCached(a:Complex,b:Complex):
#     c = getCachedComplex()
    global cacheAvail,cacheCount
    c = cacheAvail
    cacheAvail=cacheAvail.next
    cacheCount-=1
    c.r.val = a.r.val+b.r.val
    c.i.val = a.i.val+b.i.val
    return c


        
def Find_Or_Add_Complex_table(c : Complex):
    if c==cn0:
        return cn0
    if c==cn1:
        return cn1
    if abs(c.r.val-1)<epi and abs(c.i.val)<epi:
        return cn1
    if abs(c.r.val)<epi and abs(c.i.val)<epi:
        return cn0    
    
    key_r = int(c.r.val/epi)
    key_i = int(c.i.val/epi)
    res = Complex()
    if not key_r in complex_table:
        temp_r = ComplexTableEntry(c.r.val)
        complex_table[key_r] = temp_r
        res.r = temp_r
    else:
        res.r=complex_table[key_r]
    if not key_i in complex_table:
        temp_i = ComplexTableEntry(c.i.val)
        complex_table[key_i] = temp_i
        res.i = temp_i
    else:
        res.i=complex_table[key_i]
    return res

              
def getCachedComplex():
    global cacheAvail,cacheCount
    c = cacheAvail
    cacheAvail=cacheAvail.next
    cacheCount-=1
    return c

def getCachedComplex2(r_val,i_val):
    global cacheAvail,cacheCount
    c = cacheAvail
    cacheAvail=cacheAvail.next
    c.r.val = r_val
    c.i.val = i_val
    cacheCount-=1
    return c
def releaseCached(c):
    global cacheAvail,cacheCount
    c.next = cacheAvail
    cacheAvail=c
    cacheCount+=1
    return c


def equalsZero(c):
    return abs(c.r.val)<epi and abs(c.i.val)<epi

def equalsOne(c):
    return abs(c.r.val-1)<epi and abs(c.i.val)<epi

def ini_complex(max_rank=100):
    global cn0,cn1,complex_table,cacheCount,cacheAvail
    complex_table = dict()
    cacheCount = max_rank
    cacheAvail = cache_head = Complex()
#     print(cacheAvail)
    for k in range(max_rank-1):
        temp = Complex()
        cache_head.next=temp
        cache_head=temp
    
    cn0 = Find_Or_Add_Complex_table(Complex(0))
    cn1 = Find_Or_Add_Complex_table(Complex(1))
