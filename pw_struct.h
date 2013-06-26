/**
*	pw_struct.h
*	About based storage structure of  power system data, of course, some can do other things
*	Copyright(c) Peng Yao
*	Date: 2012.12
*	All rights reserved.
*/
#pragma once

#ifndef _PW_STRUCT_H
#define _PW_STRUCT_H
#endif 

#include "pw_common.h"

namespace PowerSystem
{
	template<class T>
	struct ThreeElements
	{
		int i, j;
		T data;
		ThreeElements(){}
		ThreeElements(int _i, int _j, const T& _data)
		{
			i = _i;
			j = _j;
			data = _data;
		}
	};

	struct Pair
	{
		int i;
		int j;
		Pair(int _i, int _j)
		{
			i = _i;
			j = _j; 
		}
	};

	//Undirected graph
	template<class T>
	class Graph
	{
		
	public:
		Graph();
		~Graph();
		void FromTEs(vector<ThreeElements<T>>& stEles);
		int GetNodeNum() const;
		int GetEdgeNum() const;
		int GetNodeEdgeNumById(int id) const;
		int GetNodeEdgeNumByData(const T& data) const;
		int GetNodeIdByData(const T& data) const;
		const T& GetNodeDataById(int id) const;
		void GetNodeNextById(int id, list<int>& l) const;		
		void SetNodeDataById(int id, const T& data);
		void SetConnection(int i, int j);
		void DelConnection(int i, int j);
		bool IsConnection(int i, int j) const;
		bool IdExist(int id) const;
		void AddNode(int id, const T& data);
		void DelNode(int id);
		void Clear();
		void Print() const;

		void Optimization(vector<int>& ps, int & addEdge, vector<Pair>* psAdd = NULL);

	private:
		struct node;
		typename list<node>::iterator GetNodePtr(int id);
		void GetNodeNextByPtr(typename list<node>::iterator it, list<int>& l) const;
		int GetNodeEdgeNumByPtr(typename list<node>::iterator it) const;
		void DelNodeByPtr(typename list<node>::iterator it);
	protected:
		struct node
		{
			T data;
			int id;
			list<int> first;
			node(int _id, const T& _data)
			{
				id = _id;
				data = _data;
			}
		};
		list<node> AdjList;
	};

	template<class T>
	Graph<T>::Graph()
	{
	}

	template<class T>
	void Graph<T>::Clear()
	{
		for (list<node>::iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			it->first.clear();
		}
		AdjList.clear();
	}
	
	template<class T>
	Graph<T>::~Graph()
	{
		Clear();		
	}
	
	template<class T>
	void Graph<T>::FromTEs(vector<ThreeElements<T>>& stEles)
	{
		list<Pair> ps;
		for (vector<ThreeElements<T>>::iterator it = stEles.begin(); it != stEles.end(); it++)
		{
			int i = it->i, j = it->j;
			if (i == j)
			{
				AddNode(i, it->data);
			}else if(i > j){
				Pair p(i, j);
				ps.push_back(p);

			}			
		}

		for (list<Pair>::iterator it = ps.begin(); it != ps.end(); it++)
		{
			SetConnection(it->i, it->j);
		}
	}

	template<class T>
	int Graph<T>::GetNodeNum() const
	{
		return AdjList.size();
	}

	template<class T>
	int Graph<T>::GetEdgeNum() const
	{
		int sum = 0;
		for (list<node>::iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			sum += it->first.size();
		}
		return sum;
 	}

	template<class T>
	int Graph<T>::GetNodeEdgeNumById(int id) const
	{
		if (id >= AdjList.size())
			throw _Exception("Index out of Graph's range(In GetNodeEdgeNumById)");

		for (list<node>::const_iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			if (it->id == id)
			{
				return it->first.size();
			}
		}

		return -1;
	}

	template<class T>
	int Graph<T>::GetNodeEdgeNumByPtr(typename list<node>::iterator it) const
	{
		return it->first.size();
	}

	template<class T>
	int Graph<T>::GetNodeEdgeNumByData(const T& data) const
	{
		return GetNodeEdgeNumById(GetNodeIdByData(data));
	}

	template<class T>
	int Graph<T>::GetNodeIdByData(const T& data) const
	{
		for (list<node>::const_iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			if (it->data == data)
			{
				return it->id;
			}
		}

		return T();
	}

	template<class T>
	const T& Graph<T>::GetNodeDataById(int id) const
	{
		if (id >= AdjList.size())
			throw _Exception("Index out of Graph's range(In GetNodeDataById)");

		for (list<node>::const_iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			if (it->id == id)
			{
				return it->data;
			}
		}

		return T();
	}

	template<class T>
	void Graph<T>::GetNodeNextById(int id, list<int>& l) const
	{
		GetNodeNextByPtr(GetNodePtr(id), l);	
	}

	template<class T>
	void Graph<T>::GetNodeNextByPtr(typename list<node>::iterator it, list<int>& l) const
	{
		if (it == AdjList.end())
		{
			throw _Exception("Node is not exsits(In GetNodeNextByPtr)");
		}
		l.clear();
		l.assign(it->first.begin(), it->first.end());
	}
	
	template<class T>
	void Graph<T>::SetConnection(int i, int j)
	{
		if (IsConnection(i, j))
			return;
		
		
		list<node>::iterator h = GetNodePtr(i), t = GetNodePtr(j);
		if (h == AdjList.end() || t == AdjList.end())
		{
			throw _Exception("Index out of Graph's range(In SetConnection)");
		}

		
		h->first.push_front(t->id);
		t->first.push_front(h->id);
	}

	template<class T>
	void Graph<T>::DelConnection(int i, int j)
	{
		if (!IsConnection(i, j))
			return;
	
		list<node>::iterator h = GetNodePtr(i), t = GetNodePtr(j);

		if (h == AdjList.end() || t == AdjList.end())
		{
			throw _Exception("Index out of Graph's range(In DelConnection)");
		}

		for (list<int>::iterator it = h->first.begin(); it != h->first.end();)
		{
			if (*it == t->id)
			{			
				it = h->first.erase(it);
			}else{
				it++;
			}
		}

		for (list<int>::iterator it = t->first.begin(); it != t->first.end();)
		{
			if (*it == h->id)
			{
				it = t->first.erase(it);
			}else{
				it++;
			}
		}
	}

	template<class T>
	bool Graph<T>::IsConnection(int i, int j) const
	{
		int sign = 0;
		list<node>::const_iterator h, t;
		for (list<node>::const_iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			if (it->id == i)
			{
				h = it;
				sign++;
			}
			
			if (it->id == j)
			{
				t = it;
				sign++;
			}

		}

		if (sign != 2)
			throw _Exception("Index out of Graph's range(In IsConnection)");

		for (list<int>::const_iterator it = h->first.begin(); it != h->first.end(); it++)
		{
			if (*it == j)
			{
				return true;
			}
		}
		return false;

	}

	template<class T>
	typename list<typename Graph<T>::node>::iterator Graph<T>::GetNodePtr(int id)
	{
		list<node>::iterator it;
		for (it = AdjList.begin(); it != AdjList.end(); it++)
		{
			if (it->id == id)
			{
				return it;
			}
		}
		return it;

	}
	
	template<class T>
	bool Graph<T>::IdExist(int id) const
	{
		for (list<node>::const_iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			if (it->id == id)
			{
				return true;
			}
		}
		return false;
	}

	template<class T>
	void Graph<T>::AddNode(int id, const T& data)
	{
		if (IdExist(id))
			throw _Exception("Node ID exsits(In AddNode)");
	 
		AdjList.push_back(node(id, data));
	}

	template<class T>
	void Graph<T>::DelNode(int id)
	{
		DelNodeByPtr(GetNodePtr(id));
		
	}

	template<class T>
	void Graph<T>::DelNodeByPtr(typename list<node>::iterator it)
	{
		if (it == AdjList.end())
		{
			throw _Exception("Node is not exsits(In DelNode)");
		}

		list<int> l;
		GetNodeNextByPtr(it, l);
		for (list<int>::iterator p = l.begin(); p != l.end(); p++)
		{
			DelConnection(it->id, *p);
		}
		AdjList.erase(it);
	}

	template<class T>
	void Graph<T>::SetNodeDataById(int id, const T& data)
	{
		for (list<node>::iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			if (it->id == id)
			{
				it->data = data;
			}
		}
	}

	template<class T>
	void Graph<T>::Print() const
	{
		for (list<node>::const_iterator it = AdjList.begin(); it != AdjList.end(); it++)
		{
			cout << '[' << it->id << ':' << it->data << ']';
			for (list<int>::const_iterator itj = it->first.begin(); itj!=it->first.end(); itj++)
			{
				cout << "->" << *itj;
			}
			cout << endl;
		}

	}

	//Please pay enough attention to the parameter 'addEdge' that it is the number of edges, then 1 edge means 2 injection sites
	template<class T>
	void Graph<T>::Optimization(vector<int>& ps, int & addEdge, vector<Pair>* psAdd)
	{
		ps.clear();
		
		if (psAdd != NULL)
		{
			psAdd->clear();
		}

		if (AdjList.size() < 2)
		{
			throw _Exception("Too few nodes!(In Optimization)");
		}
		ps.reserve(AdjList.size());
		addEdge = 0;
		while(AdjList.size()> 2)
		{
			list<node>::iterator it = AdjList.begin();
			int min = GetNodeEdgeNumByPtr(it);
			list<node>::iterator si = it;
			for ( ;it != AdjList.end(); it++)
			{
				if (min > GetNodeEdgeNumByPtr(it))
				{
					min = GetNodeEdgeNumByPtr(it);
					si = it;
				}
			}
			

			for (list<int>::iterator p = si->first.begin(); distance(p, si->first.end()) > 1; p++)
			{
				list<int>::iterator t = p; 
				advance(t, 1);
				for (; t != si->first.end(); t++)
				{
					
					if (!IsConnection(*p, *t))
					{
						SetConnection(*p, *t);
						addEdge++;

						if (psAdd != NULL)
						{
							Pair pair(*p, *t);
							psAdd->push_back(pair);
						}
					}
				}
			}
			ps.push_back(si->id);
			DelNodeByPtr(si);
		}
		list<node>::iterator si;
		si = AdjList.begin();
		ps.push_back(si->id);
		si++;
		ps.push_back(si->id);
	}

	template<typename T>
	class Matrix
	{
	public:
		Matrix();
		Matrix(int row);
		Matrix(int row,int col);
		Matrix(const Matrix<T> &m);
		Matrix(T** m, int row, int col, bool Static = true);
		Matrix(T* m, int row, int col);
		Matrix(T* m, int col);
		virtual ~Matrix();

		void SetSize(int row);
		void SetSize(int row,int col);		 
		void SetUnit();
		void SetValue(const T& d);
		int GetRow() const;
		int GetCol() const;

		Matrix<T> Transpose() const;
		Matrix<T> operator+(const Matrix<T>& rhs) const;
		Matrix<T> operator-(const Matrix<T>& rhs) const;
		Matrix<T> operator*(const Matrix<T>& rhs) const;
		Matrix<T> operator/(const Matrix<T>& rhs) const;
		Matrix<T> operator*(const T& d) const;
		Matrix<T> operator/(const T& d) const;
		Matrix<T>& operator=(const Matrix<T>& m);
		Matrix<T> Inverse() const;
		bool IsSquare() const; 
		bool IsSymmetric() const;
		T ToValue() const;
		T& At(int i, int j);
		T At(int i, int j) const;
		virtual void Print() const;
		void ToFile(string file = "matrixdebug.log") const;

	protected:
		void Release();
		Matrix<T> Exchange(int i,int j);
		Matrix<T> Multiple(int index, const T& mul);                  
		Matrix<T> MultipleAdd(int index,int src, const T& mul);
		int Pivot(int row) const;
		
	private:
		T** m_data;
		int Row, Col;

	};


	template<typename T>
	Matrix<T>::Matrix()
	{
		Row = 0;
		Col = 0;	
		m_data = NULL;
	}

	template<typename T>
	Matrix<T>::Matrix(int row)
	{
		Row = row;
		Col = row;

		m_data = new T*[row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
		}

	}

	template<typename T>
	Matrix<T>::Matrix(int row,int col)
	{
		Row = row;
		Col = col;

		m_data = new T*[row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
		}
	}

	template<typename T>
	Matrix<T>::Matrix(const Matrix<T> &m)
	{
		Row = m.Row;
		Col = m.Col;

		m_data = new T*[Row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
			for(int j=0;j<Col;j++)
			{
				m_data[i][j] = m.m_data[i][j];
			}

		}
	}

	template<typename T>
	Matrix<T>::Matrix(T* m, int row, int col)
	{
		Row = row;
		Col = col;

		m_data = new T*[Row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
			for(int j=0;j<Col;j++)
			{
				m_data[i][j] = *((T*)m+i*Col+j);
			}

		}
	}

	template<typename T>
	Matrix<T>::Matrix(T** m, int row, int col, bool Static)
	{
		Row = row;
		Col = col;

		m_data = new T*[Row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
			for(int j=0;j<Col;j++)
			{
				if (Static)
					m_data[i][j] = *((T*)m+i*Col+j);
				else
					m_data[i][j] = m[i][j];

			}

		}
	}

	template<typename T>
	Matrix<T>::Matrix(T* m, int col)
	{
		Row = 1;
		Col = col;

		m_data = new T*[Row];	
		m_data[0] = new T[Col];
		for(int j=0;j<Col;j++)
		{
			m_data[0][j] = m[j];

		}
	}

	template<typename T>
	Matrix<T>::~Matrix()
	{
		Release();
	}

	template<typename T>
	void Matrix<T>::Release()
	{
		if (m_data != NULL)
		{
			for (int i=0; i<Row; i++)
				delete [] m_data[i];
			delete [] m_data;
			m_data = NULL;
		}
	}

	template<typename T>
	void Matrix<T>::SetSize(int row)
	{
		Release();		
		Row = row;
		Col = row;

		m_data = new T*[row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
		}
	}

	template<typename T>
	void Matrix<T>::SetSize(int row, int col)
	{
		Release();		
		Row = row;
		Col = col;

		m_data = new T*[row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
		}
	}

	template<typename T>		 
	void Matrix<T>::SetUnit()
	{
		for(int i=0;i<Row;i++)
			for(int j=0;j<Col;j++)
				m_data[i][j] = ((i==j)?1:0);
	}

	template<typename T>
	void Matrix<T>::SetValue(const T& d)
	{
		for(int i=0;i<Row;i++)
			for(int j=0;j<Col;j++)
				m_data[i][j] = d;
	}

	template<typename T>
	int Matrix<T>::GetRow() const
	{
		return Row;
	}

	template<typename T>
	int Matrix<T>::GetCol() const
	{
		return Col;
	}

	template<typename T>
	Matrix<T> Matrix<T>::Exchange(int i,int j)
	{
		T temp;

		for(int k=0;k<Col;k++)
		{
			temp = m_data[i][k];
			m_data[i][k] = m_data[j][k];
			m_data[j][k] = temp;
		}
		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::Multiple(int index, const T& mul)    
	{
		for(int j=0;j<Col;j++)
		{
			m_data[index][j] *= mul;
		}
		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::MultipleAdd(int index, int src, const T& mul)
	{
		for(int j=0;j<Col;j++)
		{
			m_data[index][j] += m_data[src][j]*mul;
		}

		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::Transpose() const
	{
		Matrix<T> ret(Col,Row);

		for(int i=0;i<Row;i++)
			for(int j=0;j<Col;j++)
			{
				ret.m_data[j][i] = m_data[i][j];
			}
			return ret;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs) const
	{
		assert(Row == rhs.Row && Col == rhs.Col);

		Matrix<T> ret(Row, Col);

		for(int i=0;i<Row;i++)
			for(int j=0;j<Col;j++)
			{
				ret[i][j] = m_data[i][j] + rhs.m_data[i][j];
			}
			return ret;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) const
	{
		assert(Row == rhs.Row && Col == rhs.Col);

		Matrix<T> ret(Row, Col);

		for(int i=0;i<Row;i++)
			for(int j=0;j<Col;j++)
			{
				ret.m_data[i][j] = m_data[i][j] - rhs.m_data[i][j];
			}
			return ret;

	}

	template<typename T> 
	Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) const
	{
		if (Col != rhs.Row)
		{
			throw _Exception("Col != Row(In operator*)");
		}

		Matrix<T> ret(Row, rhs.Col);
		T temp;
		for(int i=0;i<Row;i++)
		{
			for(int j=0;j<rhs.Col;j++)
			{
				temp = 0;
				for(int k=0;k<Col;k++)
				{
					temp += m_data[i][k] * rhs.m_data[k][j];
				}
				ret.m_data[i][j] = temp;
			}
		}

		return ret;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator/(const Matrix<T>& rhs) const
	{
		return *this * rhs.Inverse();
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator*(const T& d) const
	{
		Matrix<T> ret(*this);
		for(int i=0;i<ret.Row;i++)
			for(int j=0;j<ret.Col;j++)
				ret.m_data[i][j] *= d;
		return ret;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator/(const T& d) const
	{
		return d*this->Inverse();
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
	{
		Release();
		Row = m.Row;
		Col = m.Col;

		m_data = new T*[Row];
		for (int i=0; i<Row; i++)
		{
			m_data[i] = new T[Col];
			for(int j=0;j<Col;j++)
			{
				m_data[i][j] = m.m_data[i][j];
			}

		}
		return *this;
	}

	template<typename T>
	int Matrix<T>::Pivot(int row) const
	{
		int index=row;

		for(int i=row+1;i<Row;i++)
		{
			if(abs(m_data[i][row]) > abs(m_data[index][row]))
				index=i;
		}

		return index;
	}

	template<typename T>
	Matrix<T> Matrix<T>::Inverse() const
	{
		if (Row != Col)
		{
			throw _Exception("Col != Row(In Inverse)");
		}

		Matrix<T> tmp(*this);
		Matrix<T> ret(Row);    
		ret.SetUnit();

		int maxIndex;
		T dMul;

		for(int i=0;i<Row;i++)
		{
			maxIndex = tmp.Pivot(i);
			if (tmp.m_data[maxIndex][i] == T(0) )
			{
				throw _Exception("0 error(In Inverse)");
			}


			if(maxIndex != i)  
			{
				tmp.Exchange(i,maxIndex);
				ret.Exchange(i,maxIndex);
			}

			ret.Multiple(i,T(1)/tmp.m_data[i][i]);
			tmp.Multiple(i,T(1)/tmp.m_data[i][i]);

			for(int j=i+1;j<Row;j++)
			{
				dMul = -tmp.m_data[j][i]/tmp.m_data[i][i];
				tmp.MultipleAdd(j,i,dMul);
				ret.MultipleAdd(j,i,dMul);

			}

		}

		for(int i=Row-1;i>0;i--)
		{
			for(int j=i-1;j>=0;j--)
			{
				dMul = -tmp.m_data[j][i]/tmp.m_data[i][i];
				tmp.MultipleAdd(j,i,dMul);
				ret.MultipleAdd(j,i,dMul);
			}
		}       
		return ret;
	}

	template<typename T>
	bool Matrix<T>::IsSquare() const
	{
		return Row==Col;
	}

	template<typename T>
	bool Matrix<T>::IsSymmetric() const
	{

		if(Row != Col) 
			return false;

		for(int i=0;i<Row;i++)
			for(int j=i+1;j<Col;j++)
				if( m_data[i][j] != m_data[j][i])
					return false;

		return true;
	}

	template<typename T>
	T Matrix<T>::ToValue() const
	{
		if (Row != 1 || Col != 1)
		{
			throw _Exception("Row != 1 or Col != 1(In ToValue)");
		}
		return m_data[0][0];
	}

	template<typename T>
	T& Matrix<T>::At(int i, int j)
	{
		if (i >= Row || j >= Col)
			throw _Exception("i >= Row or j >= Col(In At)");
		return m_data[i][j];
	}

	template<typename T>
	T Matrix<T>::At(int i, int j) const
	{
		if (i >= Row || j >= Col)
			throw _Exception("i >= Row or j >= Col(In At)");
		return m_data[i][j];
	}
	template<typename T>
	void Matrix<T>::Print() const
	{	
		for(int i=0;i<Row;i++)
		{
			for(int j=0;j<Col;j++)
			{
				cout << m_data[i][j] << ',';
			}
 			cout << "\r\n";
		}
		cout << "\r\n";
	 
	}

	template<typename T>
	void Matrix<T>::ToFile(string file) const
	{
		fstream fp;
		fp.open(file.c_str(), ios::out | ios::app);
		for(int i=0;i<Row;i++)
		{
			for(int j=0;j<Col;j++)
			{
				fp << m_data[i][j] << ',';
			}
 			fp << "\r\n";
		}
		fp << "\r\n";
		fp.flush();
		fp.clear();
		fp.close();
	}


	template<class T>
	class SparseMatrix
	{
	public:
		SparseMatrix();
		SparseMatrix(vector<ThreeElements<T>>& stEles, int total, int rows);
		SparseMatrix(const SparseMatrix<T>& sm);
		~SparseMatrix();
		SparseMatrix<T> & operator=(const SparseMatrix<T> & sm);
	
		void DelElement(int i, int j);
		bool IsZero(int i, int j);
		T GetAt(int i, int j) const;
		void GetAt(int i, int j, T* &data) const;
		void SetAt(int i, int j, const T& data);
		void SetAt(const ThreeElements<T>& stEle);

		virtual void InitFromFile(string strfile);
		virtual void ToGraph(Graph<T>& gh) const;
		virtual void Clear();
		virtual string ToString() const;
		void ToMatrix(Matrix<T>& matrix) const;
		void ToTEs(vector<ThreeElements<T>>& stEles, int &rows) const;
		void FromMatrix(const Matrix<T>& matrix);
		void FromTEs(vector<ThreeElements<T>>& stEles, int total, int rows);
		void Print() const;
		void LDUDecomposition(bool noInjection = true, int layer = 0);
		void OptimizeLDU();
		static void SortTEs(vector<ThreeElements<T>>& stEles, int rc);
		static void MatrixToTEs(const Matrix<T>& matrix, vector<ThreeElements<T>>& stEles, int &rows);
		void SolveEquation(vector<T> & sm);

	protected:
		void AddElement(int i, int j, const T& data);
		void AddElement(const ThreeElements<T>& stEle);
		int LINKBeforeAt(int i, int j, bool * isFirst = NULL) const;
		bool SumMulLdu(int i, int j, int k, T& ans) const;
		
		
	protected:
		vector<T> VA;			//store non-zero elements in the matrix by row.
	private:
		vector<int> JA;			//store column number of non-zero elements in the matrix by row. 
		vector<int> LINK;		//store the position of next non-zero element in VA, the last non-zero element of each line is replaced by 0.
		vector<int> IA;			//record the position of first non-zero element of each line in VA.
		vector<int> NA;			//record the number of non-zero elements of each line.
	};

	template<class T>
	SparseMatrix<T>::SparseMatrix()
	{
	}

	template<class T>
	SparseMatrix<T>::SparseMatrix(vector<ThreeElements<T>>& stEles, int total, int rows)
	{		
		VA.reserve(total);
		JA.reserve(total);
		
		LINK.resize(total, -1);	
		IA.resize(rows);
		NA.resize(rows);
		
		//SortTEs(stEles, rows);  //Call it by hand is better
		for (vector<ThreeElements<T>>::const_iterator it = stEles.begin(); it != stEles.end(); it++)
		{
			int i = it->i, j = it->j;
			T data = it->data;
			VA.push_back(data);
			JA.push_back(j);
			NA[i]++;
		}
	
		for (int i=0; i<rows; i++)
		{
			int sum = 0;
			for (int k=0; k<i; k++)
			{
				sum += NA[k];
			}
			IA[i] = sum;
			for (int j=0; j<NA[i] - 1; j++)
			{
				LINK[sum + j] = sum + j + 1;

			}
		}
	}

	template<class T>
	void SparseMatrix<T>::FromTEs(vector<ThreeElements<T>>& stEles, int total, int rows)
	{	
		Clear();
		VA.reserve(total);
		JA.reserve(total);
		
		LINK.resize(total, -1);	
		IA.resize(rows);
		NA.resize(rows);
		
		for (vector<ThreeElements<T>>::const_iterator it = stEles.begin(); it != stEles.end(); it++)
		{
			int i = it->i, j = it->j;
			T data = it->data;
			VA.push_back(data);
			JA.push_back(j);
			NA[i]++;
		}
	
		for (int i=0; i<rows; i++)
		{
			int sum = 0;
			for (int k=0; k<i; k++)
			{
				sum += NA[k];
			}
			IA[i] = sum;
			for (int j=0; j<NA[i] - 1; j++)
			{
				LINK[sum + j] = sum + j + 1;

			}
		}
	}
	
	template<class T>
	SparseMatrix<T>::SparseMatrix(const SparseMatrix& sm)
	{
		VA =sm.VA;
		JA = sm.JA;
		LINK = sm.LINK;
		IA = sm.IA;
		NA = sm.NA;
	}
	
	template<class T>
	SparseMatrix<T>::~SparseMatrix()
	{
		Clear();
	}

	template<class T>
	SparseMatrix<T> & SparseMatrix<T>::operator=(const SparseMatrix<T> & sm)
	{
		Clear();
		VA.assign(sm.VA.begin(), sm.VA.end());
		JA.assign(sm.JA.begin(), sm.JA.end());
		LINK.assign(sm.LINK.begin(), sm.LINK.end());
		IA.assign(sm.IA.begin(), sm.IA.end());
		NA.assign(sm.NA.begin(), sm.NA.end());
		return *this;
	}

	template<class T>
	void SparseMatrix<T>::ToGraph(Graph<T>& gh) const
	{
		int rows = NA.size();

		list<Pair> ps;
		for (int i=0; i<rows; i++)
		{
			int k = IA[i];
			while (k != -1)
			{
				int j = JA[k];

				if (i == j)
				{
					gh.AddNode(i, VA[k]);
				}else if(i > j){
					Pair p(i, j);
					ps.push_back(p);

				}
				k = LINK[k];
			}
		}

		for (list<Pair>::iterator it = ps.begin(); it != ps.end(); it++)
		{
			gh.SetConnection(it->i, it->j);
		}

	}

	template<class T>
	void SparseMatrix<T>::Clear()
	{
		VA.clear();
		JA.clear();
		LINK.clear();
		IA.clear();
		NA.clear();
	}

	template<class T>
	void SparseMatrix<T>::AddElement(int i, int j, const T& data)
	{
		VA.push_back(data);
		JA.push_back(j);
		bool first;
		int k = LINKBeforeAt(i, j, &first);
		LINK.push_back(LINK[k]);
		int last = VA.size() - 1;
		if (first)
		{
			IA[i] = last;
			LINK[last] = k;
		}else{
			LINK[last] = LINK[k];
			LINK[k] = last;		
		}
		
		NA[i]++;
	}

	template<class T>
	void SparseMatrix<T>::AddElement(const ThreeElements<T>& stEle)
	{
		AddElement(stEle.i, stEle.j, stEle.data);
	}

	template<class T>
	int SparseMatrix<T>::LINKBeforeAt(int i, int j, bool* isFirst) const
	{
		if ((unsigned)i >= IA.size())
			throw _Exception("Index is out of range(In LINKBeforeAt)");
		int k = IA[i];
		int temp = k;

		bool first = true;
		while (k != -1)
		{
			int _j = JA[k];
			if (j <= _j)
				break;
			first = false;
			temp = k;
			k = LINK[k];
			
		}

		if (isFirst != NULL)
			*isFirst = first;
		return temp;
	}

	template<class T>
	void SparseMatrix<T>::DelElement(int i, int j)
	{
		if (IsZero(i, j))
		{
			return;
		}

		bool first;
		int k = LINKBeforeAt(i, j, &first);
		if (first)
		{
			IA[i] = LINK[k];
		}else{
			LINK[k] = LINK[LINK[k]];
		}
		
		NA[i]--;
		
	}

	template<class T>
	bool SparseMatrix<T>::IsZero(int i, int j)
	{
		int k = IA[i];
		
		while (k != -1)
		{
			int _j = JA[k];

			if (j == _j)
				return false;

			if (j < _j)
				break;
			k = LINK[k];
		}

		return true;
	}

	template<class T>
	T SparseMatrix<T>::GetAt(int i, int j) const
	{
		int k = IA[i];

		while (k != -1)
		{
			int _j = JA[k];

			if (j == _j)
				return VA[k];

			if (j < _j)
				break;
			k = LINK[k];
		}
		return T(0.0);
	}

	template<class T>
	void SparseMatrix<T>::GetAt(int i, int j, T* &data) const
	{
		int k = IA[i];

		while (k != -1)
		{
			int _j = JA[k];

			if (j == _j)
			{
				data = &VA[k];
				return;
			}

			if (j < _j)
				break;
			k = LINK[k];
		}
		data = NULL;
	}

	template<class T>
	void SparseMatrix<T>::SetAt(int i, int j, const T& data)
	{
		if (IsZero(i, j))
		{
			AddElement(i, j, data);
		}else{
			int k = IA[i];

			while (k != -1)
			{
				int _j = JA[k];

				if (j == _j)
				{
					VA[k] = data;
					break;
				}
				k = LINK[k];
					
			}
		}
	}

	template<class T>
	void SparseMatrix<T>::SetAt(const ThreeElements<T>& stEle)
	{
		SetAt(stEle.i, stEle.j, stEle.data);
	}

	template<class T>
	void SparseMatrix<T>::InitFromFile(string strfile)
	{
		fstream fp;
		fp.open(strfile.c_str(), ios::in);

		if (!fp.is_open())
			throw _Exception(string("Can not open file[") + strfile + string("](In InitFromFile)"));

		fp.clear();
		fp.close();

	}

	template<class T>
	void SparseMatrix<T>::SortTEs(vector<ThreeElements<T>>& stEles, int rc)
	{
		for (vector<ThreeElements<T>>::iterator it = stEles.begin(); it != stEles.end() - 1; it++)
		{
			for(vector<ThreeElements<T>>::iterator itj = it + 1; itj != stEles.end(); itj++)
			{
				if (it->i * rc + it->j > itj->i * rc + itj->j )
				{
					ThreeElements<T> temp;
					temp = *it;
					*it = *itj;
					*itj = temp;
				}
			}
		}	
	}

	template<class T>
	string SparseMatrix<T>::ToString() const
	{
		return string("");
	}

	template<class T>
	void SparseMatrix<T>::ToMatrix(Matrix<T> & matrix) const
	{
		int rows = NA.size();

		matrix.SetSize(rows);
		matrix.SetValue(T(0.0));

		for (int i=0; i<rows; i++)
		{
			int k = IA[i];
			while (k != -1)
			{
				int j = JA[k];
				matrix.At(i, j) = VA[k];
				k = LINK[k];
			}
		}
	}
	
	template<class T>
	void SparseMatrix<T>::FromMatrix(const Matrix<T>& matrix)
	{
		vector<ThreeElements<T>> stEles; 
		int rows;
		MatrixToTEs(matrix, stEles, rows);
		
		//No Efficiency!
		int total = stEles.size();
		VA.reserve(total);
		JA.reserve(total);
		
		LINK.resize(total, -1);
		IA.resize(rows);
		NA.resize(rows);
	
		
		for (vector<ThreeElements<T>>::const_iterator it = stEles.begin(); it != stEles.end(); it++)
		{
			int i = it->i, j = it->j;
			T data = it->data;
			VA.push_back(data);
			JA.push_back(j);
			NA[i]++;
		}
	
		for (int i=0; i<rows; i++)
		{
			int sum = 0;
			for (int k=0; k<i; k++)
			{
				sum += NA[k];
			}
			IA[i] = sum;
			for (int j=0; j<NA[i] - 1; j++)
			{
				LINK[sum + j] = sum + j + 1;

			}
		}
		
	}

	template<class T>
	void SparseMatrix<T>::MatrixToTEs(const Matrix<T>& matrix, vector<ThreeElements<T>>& stEles, int &rows)
	{
		if (!matrix.IsSquare())
			throw _Exception("Not squre matrix(In MatrixToTEs)");
		
		rows = matrix.GetRow();
		stEles.clear();
		for (int i=0; i<matrix.GetRow(); i++)
			for (int j=0; j<matrix.GetCol(); j++)
			{
				if (fabs(matrix.At(i, j) - T(0.0)) <= T(_ZERO))
					continue;
				ThreeElements<T> temp(i, j, matrix.At(i, j));
				stEles.push_back(temp);
			}
	}

	template<class T>
	void SparseMatrix<T>::ToTEs(vector<ThreeElements<T>>& stEles, int &rows) const
	{
		rows = NA.size();
		stEles.clear();
		for (int i=0; i<rows; i++)
		{
			int k = IA[i];
			while (k != -1)
			{
				int j = JA[k];

				ThreeElements<T> temp(i, j, VA[k]);
				stEles.push_back(temp);
				k = LINK[k];
			}
		}
	}

	template<class T>
	void SparseMatrix<T>::Print() const
	{
		vector<T>::const_iterator ciT;
		vector<int>::const_iterator ciI;
		cout << "VA:";
		for (ciT = VA.begin(); ciT != VA.end(); ciT++)
		{
			cout << *ciT << "   ";
		}
		cout << endl;

		cout << "JA:";
		for (ciI = JA.begin(); ciI != JA.end(); ciI++)
		{
			cout << *ciI << "   ";
		}
		cout << endl;

		cout << "LK:";
		for (ciI = LINK.begin(); ciI != LINK.end(); ciI++)
		{
			cout << *ciI << "   ";
		}
		cout << endl;

		cout << "IA:";
		for (ciI = IA.begin(); ciI != IA.end(); ciI++)
		{
			cout << *ciI << "   ";
		}
		cout << endl;

		cout << "NA:";
		for (ciI = NA.begin(); ciI != NA.end(); ciI++)
		{
			cout << *ciI << "   ";
		}
		cout << endl;
	}

	/*
	Like this:
				  _j
	
		d---------u
		|  d------u
		|  |  d---u
		|  |  |
	_i	l  l  l  +(_k)
	*/
	template<class T>
	bool SparseMatrix<T>::SumMulLdu(int _i, int _j, int _k, T& ans) const
	{
		ans = T(0.0);
		//ans = 0 -> false
		bool ret = false;
		int k = IA[_i];
		while (k != -1)
		{
			//Find li
			int j = JA[k];
			T l = VA[k], d, u;
			if (j > _k)
			{
				break;  
			}
			//Find uj
			int kk = IA[j];
			while (kk != -1)
			{
				int jj = JA[kk];
				// find d
				if (jj == j)
				{
					d = VA[kk];
				}
				if (jj == _j)
				{
					u = VA[kk];
					ans += l * d * u;
					ret = true;
					break;
				}
				if (jj > _j)
					break;
				kk = LINK[kk];
			}

			k = LINK[k];
		}
		return ret;

	}

	//Assuming strictly diagonally dominant matrix
	template<class T>
	void SparseMatrix<T>::LDUDecomposition(bool noInjection, int layer)
	{
		if ((unsigned)layer >= NA.size())
		{
			OptimizeLDU();
			return;
		}
		if (0 == layer)
		{
			for (int k = 1; k < (int)JA.size(); k++)
			{
				if (JA[k] == 0)
				{
					VA[k] /= VA[0];
				}
				if (k < NA[0])
				{
					VA[k] /= VA[0];
				}
			}
			LDUDecomposition(noInjection, ++layer);
		}else{
			int k = IA[layer];
			int i = layer;
			T d, ans;
			int pursuit = i + 1;
			while (k != -1)
			{
				int j = JA[k];
				
				if (i == j)
				{
					if (SumMulLdu(i, i, i - 1, ans))
					{
						VA[k] = VA[k] -  ans;
					}
					d = VA[k];
				}else if(i < j)
				{
					if (SumMulLdu(i, j, i - 1, ans))
					{
						VA[k] = (VA[k] -  ans);
					}
					VA[k] /= d;

					if (!noInjection)
					{
						
						for (; pursuit < j; pursuit++)
						{
							if (SumMulLdu(i, pursuit, i - 1, ans))
							{
								//Inject by row
								AddElement(i, pursuit, ans / d);
							}
						}
						pursuit = j + 1;
					}
				}		
				k = LINK[k];
			}
			
			if (!noInjection)
			{			
				for (; pursuit  < (int)NA.size(); pursuit++)
				{
					if (SumMulLdu(i, pursuit, i - 1, ans))
					{
						//Inject by row
						AddElement(i, pursuit, (T(0.0) - ans) / d);
					}
				}
			}

			for (int _i = i + 1; _i < (int)NA.size(); _i++)
			{
				k = IA[_i];
				bool sign = false;
				while (k != -1)
				{
					if (JA[k] == i)
					{
						if (SumMulLdu(_i, i, i - 1, ans))
						{
							VA[k] = (VA[k] -  ans);
						}
						VA[k] /= d;
						sign = true; //No injection
						break;
					}else if(JA[k] >= i)
					{
						break;
					}
					k = LINK[k];
				}

				if (!noInjection && !sign)
				{
					if (SumMulLdu(_i, i, i - 1, ans))
					{
						//Inject by col
						AddElement(_i, i, (T(0.0) - ans) / d);
					}
				}			
			}	
			LDUDecomposition(noInjection, ++layer);
		}
	}

	template<class T>
	void SparseMatrix<T>::OptimizeLDU()
	{
		for (int i=0; i < (int)NA.size(); i++)
		{
			int k = IA[i];

			while (k != -1)
			{
				int j = JA[k];

				if (i == j)
				{
					VA[k] = T(1.0) / VA[k];
					break;
				}else if (j > i)
				{
					break;
				}
				k = LINK[k];
			}
		}
	}

	template<class T>
	void SparseMatrix<T>::SolveEquation(vector<T> & sm)
	{
		if (NA.size() != sm.size())
			throw _Exception("the size is different(In SolveEquation)");

		//forward substitution
		for (int i = 0; i < (int)NA.size() - 1; i++)
		{
			if (fabs(sm[i] - T(0.0)) > _ZERO)
			{
				for (int j = i + 1; j < (int)NA.size(); j++)
				{
					int k = IA[j];

					while (k != -1)
					{
						if (JA[k] == i)
						{
							sm[j] -= VA[k] * sm[i];
							break;
						}else if (JA[k] > i)
						{
							break;
						}
						k = LINK[k];
					}
				}
			}
		}

		//Normalized
		for (int i = 0; i < (int)NA.size(); i++)
		{
			if (fabs(sm[i] - T(0.0)) > _ZERO)
			{		
				int k = IA[i];

				while (k != -1)
				{
					if (JA[k] == i)
					{
						sm[i] *= VA[k];
						break;
					}else if (JA[k] > i)
					{
						break;
					}
					k = LINK[k];
				}
			}
		}

		//backward substitution
		for (int j = (int)NA.size() - 1; j > 0 ; j--)
		{
			if (fabs(sm[j] - T(0.0)) > _ZERO)
			{
				for (int i = j - 1; i >= 0 ; i--)
				{
					int k = IA[i];

					while (k != -1)
					{
						if (JA[k] == j)
						{
							sm[i] -= VA[k] * sm[j];
							break;
						}
						k = LINK[k];
					}
				}
			}
		}
	}
};