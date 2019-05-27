import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import nju.lzx.Algorithm.Greedy;
import nju.lzx.Algorithm.GreedyGeneral;
import nju.lzx.Algorithm.TabuSearch;
import nju.lzx.Constraint.CapacityConstraint;
import nju.lzx.Constraint.MinimizeDistance;
import nju.lzx.Constraint.MinimizeFixedCost;
import nju.lzx.Data.Instance;
import nju.lzx.InitialSolutionOperator.InsertBase;
import nju.lzx.Interface.Constraint;
import nju.lzx.Interface.Operator;
import nju.lzx.Interface.Route;
import nju.lzx.LocalSearchOperator.CrossBase;
import nju.lzx.LocalSearchOperator.ExchangeBaseDeep;
import nju.lzx.LocalSearchOperator.RelocateBase;
import nju.lzx.LocalSearchOperator.RelocateBaseIntra;
import nju.lzx.Route.RouteBase;
import nju.lzx.Utility.Atr;
import nju.lzx.Utility.Reference;


public class Main {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		double t1 = System.nanoTime();
		Instance inst = load_instance("data/E-n101-k14.vrp");
		//设置模式参数
		inst.parameter.Mode.multi_thread_enable = false;
		//设置初始解参数
		inst.parameter.InitialSolution.log_print = false;
		//设置禁忌搜索参数
		inst.parameter.TabuSearch.maximum_iteration = 2000;
		inst.parameter.TabuSearch.maximum_tabu_tenure = 50;
		inst.parameter.TabuSearch.tenure_decay_rate = 0.99;
		inst.parameter.TabuSearch.mininum_tabu_tenure = 10;
		inst.parameter.TabuSearch.mininum_shake_tenure = 100;
		inst.parameter.TabuSearch.minimum_shake_iteration = 100;
		inst.parameter.TabuSearch.log_print = true;
		inst.parameter.TabuSearch.log_detail = false;
		//设置算子参数
		inst.parameter.Operator.insertion_prune_threshhold = 1e6;
		inst.parameter.Operator.exchange_prune_threshhold = 1e6;
		inst.parameter.Operator.cross_prune_threshhold = 1e6;
		inst.parameter.Operator.remove_prune_threshhold = 0;
		inst.parameter.Operator.route_cross_threshhold = 1e6;
		
		
		//构造约束条件
		Constraint[] cnts = new Constraint[2];
		MinimizeFuelCost.ConstraintData dat = new MinimizeFuelCost.ConstraintData(inst.d, inst.q, 1.51, 17.22);
		cnts[0] = new MinimizeFuelCost(dat, 1);
		CapacityConstraint.ConstraintData cap_dat = new CapacityConstraint.ConstraintData(inst.q, inst.Q);
		cnts[1] = new CapacityConstraint(cap_dat, true, 100);
		
		//构造算法算子
		Operator[] operators = new Operator[4];
		double[] coefs = new double[4];
		operators[0] = new RelocateBase(inst);
		coefs[0] = 1;
		operators[1] = new ExchangeBaseDeep(inst);
		coefs[1] = 0.5;
		operators[2] = new CrossBase(inst);
		coefs[2] = 1;
		operators[3] = new RelocateBaseIntra(inst);
		coefs[3] = 1;
		
		//构造需要访问的节点集合
		ArrayList<Atr> atrs = new ArrayList<Atr>();
		for(int i = 1; i < inst.n; i++){
			atrs.add(new Atr(i));
		}
		boolean[] exc = new boolean[inst.n];
		for(int i = 1; i < inst.n; i++)
			exc[i] = true;
		
		//构造初始解
		Greedy greedy = new Greedy(inst, new InsertBase(inst, cnts));
		ArrayList<Route> s = greedy.generate(atrs);
		System.out.println(greedy.toString(s));
		greedy.check(s, true, exc);
		System.out.println("feasibility of the initial solution>>>" + greedy.is_feasible(s) + "\t" + s.size() + "\t" + greedy.get_total_cost(s));
		
		//最小化行驶距离
		TabuSearch tabu = new TabuSearch(inst, operators, coefs);
		s = toDeep(inst, s);
		tabu.check(s, true, exc);
		for(int i = 0; i < s.size(); i++){
			//s.get(i).relax(1);
		}
		ArrayList<Route> bs = tabu.solve(s);

		tabu.check(bs, true, exc);
		//System.out.println(tabu.toString(bs));
		double t2 = System.nanoTime();
		System.out.println("computation time>>>" + (t2 - t1) / 1e9);
		//inst.statistics.print();
	}
	

	public static Instance load_instance(String path) throws FileNotFoundException{
		Instance inst = new Instance();
		Scanner cin = new Scanner(new BufferedReader(new FileReader(path)));
		for(int i = 0; i < 3; i++)
			cin.nextLine();
		for(int i = 0; i < 2; i++)
			cin.next();
		inst.n = cin.nextInt();
		inst.q = new double[inst.n];
		inst.lng = new double[inst.n];
		inst.lat = new double[inst.n];
		inst.d = new double[inst.n][inst.n];
		inst.t = new double[inst.n][inst.n];
		for(int i = 0; i < 5; i++)
			cin.next();
		inst.Q = cin.nextDouble();
		cin.next();
		for(int i = 0; i < inst.n; i++){
			cin.next();
			inst.lng[i] = cin.nextDouble();
			inst.lat[i] = cin.nextDouble();
		}
		cin.next();
		for(int i = 0; i < inst.n; i++){
			cin.next();
			inst.q[i] = cin.nextDouble();
		}
		for(int i = 0; i < inst.n; i++){
			for(int j = i + 1; j < inst.n; j++){
				inst.d[i][j] = inst.d[j][i] = inst.t[i][j] = inst.t[j][i] = 
						Math.sqrt((inst.lng[i] - inst.lng[j]) * (inst.lng[i] - inst.lng[j]) + (inst.lat[i] - inst.lat[j]) * (inst.lat[i] - inst.lat[j]));
			}
		}
		cin.close();
		return inst;
	}

	
	public static ArrayList<Route> toDeep(Instance inst, ArrayList<Route> s){
		ArrayList<Route> ns = new ArrayList<Route>();
		for(int i = 0; i < s.size(); i++){
			Route r = s.get(i);
			Reference ref = r.get_reference();
			Route nr = new RouteBase(inst, ref.len, ref.seq, false, true);
			Constraint[] _cnts = r.get_constraints();
			Constraint[] _cnts2 = new Constraint[_cnts.length];
			for(int j = 0; j < _cnts.length; j++){
				_cnts2[j] = _cnts[j].copy(nr.get_reference());
			}
			nr.add_constraints(_cnts2);
			ns.add(nr);
		}
		return ns;
	}
}
