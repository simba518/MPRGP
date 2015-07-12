#ifndef READ_MATLAB_PROB_H
#define READ_MATLAB_PROB_H

#include "MatVec.h"

PRJ_BEGIN

static void readMatlabProb(std::istream& is,Kernel<scalarD>::Vec& b,
                           FixedSparseMatrix<scalarD,Kernel<scalarD> >& A)
{
    int row,nnz;
    readBinaryData(row,is);
    readBinaryData(nnz,is);

    A.resize(row,row);
    b.resize(row);

    std::vector<int> rowStart(row+1);
    for(sizeType i=0;i<row+1;i++)
        readBinaryData(rowStart[i],is);

    std::vector<std::pair<scalarD,int> > elem(nnz);
    for(sizeType i=0;i<nnz;i++)
    {
        readBinaryData(elem[i].second,is);
        readBinaryData(elem[i].first,is);
    }

    std::vector<Eigen::Triplet<scalarD,sizeType> > triplets;
    for(sizeType i=0;i<row;i++)
    for(sizeType j=rowStart[i];j<rowStart[i+1];j++)
        triplets.push_back(Eigen::Triplet<scalarD,sizeType>(i,elem[j].second,elem[j].first));
    A.buildFromTriplets(triplets.begin(),triplets.end());
    ASSERT(A.isSymmetric(1E-9f));

    b.resize(row);
    for(sizeType i=0;i<row;i++)
        readBinaryData(b[i],is);
}

PRJ_END

#endif