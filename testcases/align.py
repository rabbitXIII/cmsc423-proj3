
Infinity = float('inf')

def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in xrange(sizex)]

def print_matrix(x, y, A):

    # print the top row
    print "%5s %5s" % (" ","*"),
    for c in x:
        print "%5s" % (c),
    print

    y = "*" + y
    for j in xrange(len(A[0])):
        print "%5s" % (y[j]),
        for i in xrange(len(A)):
            print "%5.0f" % (A[i][j]),
        print


class ScoreParam:
    """Stores the parameters for an alignment scoring function"""
    def __init__(self, gap_start, gap, match, mismatch):
        self.gap_start = gap_start
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

    def matchchar(self, a,b):
        return self.match if a==b else self.mismatch

    def __str__(self):
        return "match = %d; mismatch = %d; gap_start = %d; gap_extend = %d" % (self.match, self.mismatch, self.gap_start, self.gap)


def local_align(x, y, score=ScoreParam(0, -7, 10, -5)):
    """Do a local alignment between x and y"""
    # create a zero-filled matrix
    A = make_matrix(len(x) + 1, len(y) + 1)

    best = 0
    optloc = (0,0)

    # fill in A in the right order
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):

            # the local alignment recurrance rule:
            A[i][j] = max(
               A[i][j-1] + score.gap,
               A[i-1][j] + score.gap,
               A[i-1][j-1] + (score.match if x[i-1] == y[j-1] else score.mismatch),
               0
            )

            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j)

    print_matrix(x, y, A)
    # return the opt score and the best location
    return best, optloc


def affine_align(x, y, score=ScoreParam(-15, -7, 10, -2)):
    """Global alignment with affine penalties"""
    M = make_matrix(len(x) + 1, len(y) + 1)
    X = make_matrix(len(x) + 1, len(y) + 1)
    Y = make_matrix(len(x) + 1, len(y) + 1)

    for i in xrange(1, len(x)+1):
        M[i][0] = -Infinity
        X[i][0] = -Infinity
        Y[i][0] = score.gap_start + i * score.gap

    for i in xrange(1, len(y)+1):
        M[0][i] = -Infinity
        X[0][i] = score.gap_start + i * score.gap
        Y[0][i] = -Infinity

    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):

            M[i][j] = score.matchchar(x[i-1], y[j-1]) + max(
                    M[i-1][j-1],
                    X[i-1][j-1],
                    Y[i-1][j-1]
            )

            X[i][j] = max(
                    score.gap_start + score.gap + M[i][j-1],
                    score.gap + X[i][j-1],
                    score.gap_start + score.gap + Y[i][j-1]
            )

            Y[i][j] = max(
                    score.gap_start + score.gap + M[i-1][j],
                    score.gap_start + score.gap + X[i-1][j],
                    score.gap + Y[i-1][j]
            )

    opt = max(M[len(x)][len(y)], X[len(x)][len(y)], Y[len(x)][len(y)])

    print "x = %s & y = %s" % (x,y)
    print "Scoring:", str(score)
    print "M matrix ="
    print_matrix(x,y,M)
    print "X matrix ="
    print_matrix(x,y,X)
    print "Y matrix ="
    print_matrix(x,y,Y)
    print "Optimal =", opt

    return opt



if __name__ == "__main__":
    x = "AGCAGGGGT"
    y = "CAGG"
    #                   (gap_start, gap, match, mismatch)
    affine_align(x, y, score=ScoreParam(-15, -7, 30, -5))

