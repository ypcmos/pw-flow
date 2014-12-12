/**
*	pw_flow.h
*	About power flow calculation in power system
*	Copyright(c) Peng Yao
*	Date: 2012.12
*	All rights reserved
*/
#pragma once

#ifndef _PW_FLOW_H
#define _PW_FLOW_H
#endif 

#include "pw_struct.h"
#include <complex>

namespace PowerSystem
{
	class BusAndBranch
	{
	public:	 
		struct Branch 
		{
			int	i;   
			int j;    
			string Name;
			enum BranchType {LINE = 1, TRANSFORMER};
			BranchType Type;
			double R;
			double X;  
			double B;  
			double K; 
			Branch(){}
			Branch(int _i, int _j, string _Name, BranchType _Type, double _R, double _X, double _B, double _K)
			{
				i = _i;
				j = _j;
				Name = _Name;
				Type = _Type;
				R = _R;
				X = _X;
				B = _B;
				K = _K;
			}
		};

		struct Bus
		{
			int i;
			string Name;
			enum NodeType {SLACK = 1, PV, PQ}; 
			NodeType Type;
			double V,cita, BS, GP, GQ, LP, LQ, Qmax, Qmin, b, P, Q;
			Bus(){}
			Bus(int _i, string _Name, NodeType _Type, double _V, double _cita, double _BS, double _GP, double _GQ, double _LP, double _LQ, double _Qmax, double _Qmin, double _b)
			{ 
				i = _i;
				Name = _Name;
				Type = _Type;
				V = _V;
				cita = _cita;
				BS = _BS;
				GP = _GP / BS;
				GQ = _GQ / BS;
				LP = _LP / BS;
				LQ = _LQ / BS;
				Qmax = _Qmax / BS;
				Qmin = _Qmin / BS;
				b = _b;
				P = 0.0;
				Q = 0.0;
			}
		};
	
	public:
		BusAndBranch();
		void FromConfigFile(string filename);
		int GetBusSize() const;
		int GetBranchSize() const;
		Bus* GetBus(int id);
		Branch* GetBranch(int id);
		static void PrintBranch(const Branch& one);	
		static void PrintBus(const Bus& one);
		void PrintBranches() const;
		void PrintBuses() const;
		
	private:
		vector<Branch> branch;
		vector<Bus> bus;
	};

	class _PowerFlow
	{
	public:
		_PowerFlow();
		virtual ~_PowerFlow();
		virtual void CreateYMatrix() = 0;
		virtual void InitBABFromFile(string);	
		virtual void SortBABByType();
		virtual int GetSystemMap(int i) const;
		virtual int GetVirtualId(int i);
		virtual int GetNodeNum() const;
		virtual void FlowCal(int type = 0, int times = 10, double e = 1.0e-5) = 0;
		virtual void BranchFlowCal();
		virtual double GetPi(int i, const vector<double>& V, const vector<double>& cita) const;
		virtual double GetQi(int i, const vector<double>& V, const vector<double>& cita) const;
	protected:
		BusAndBranch BAB;
		vector<Pair> NodeIdMap; 
		vector<ThreeElements<complex<double>>> Y;
		int nSlack, nPQ, nPV;
		vector<ThreeElements<double>>Pij, Qij;
		int Times;
	};

	class FastPowerFlow:public _PowerFlow
	{
	public:
		enum Type {XB = 1, BX, PQ};
	public:
		FastPowerFlow();
		virtual ~FastPowerFlow();
		virtual void CreateYMatrix();	
		virtual void FlowCal(int type = 0, int times = 10, double e = 1.0e-5);
		virtual void CreateB1_XB();
		virtual void CreateB2_XB();
		virtual void CreateB1_BX();
		virtual void CreateB2_BX();
		virtual void CreateB1_PQ();
		virtual void CreateB2_PQ();
		virtual void GetB1FactorTable();
		virtual void GetB2FactorTable();
		virtual void CreateB1();
		virtual void CreateB2();
		virtual void Print() const;
		virtual void ToFile(string file = "flow.log");

	private:
		
		int GetB2Virtual(int i) const;
		int GetB1Virtual(int i) const;

	private:
		SparseMatrix<double>B1, B2;
		vector<int> B1Map, B2Map;
		Type type;
		
	};

	class NRPowerFlow: public _PowerFlow
	{
	public:
		NRPowerFlow();
		virtual ~NRPowerFlow();
		virtual void CreateYMatrix();	
		virtual void FlowCal(int type = 0, int times = 10, double e = 1.0e-4);
		virtual void Print() const;
		virtual void ToFile(string file = "flow.log");

	private:		
		void InitJacobi();
		int GetDistance(int i, int j);
	private:
		int JacobiSize;
		vector<Pair> JacobiPos;
		vector<int> JacobiDis;
	
	};
	void test_Print(const vector<int>& v);
	void test_Print(const vector<double>& v);
	void test_Print(const vector<Pair>& v);
	void test_Print(const vector<ThreeElements<double>>& v);
};