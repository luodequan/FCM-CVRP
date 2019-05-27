
import nju.lzx.Interface.Base;
import nju.lzx.Interface.Constraint;
import nju.lzx.Constraint.AbstractConstraint;
import nju.lzx.Utility.Reference;

public class MinimizeFuelCost extends AbstractConstraint implements Base {

	/** 约束的问题数据。 */
	protected ConstraintData dat;
	
	/** fq[i]：从仓库到第i个节点的需求之和。*/
	protected double[] fq;
	
	/** bd[i]：从第i个结点回到仓库的距离。 */
	protected double[] bd;
	
	
	public static class ConstraintData{
		
		/** 节点需求。 */
		public double[] q;
		
		/** 距离矩阵。 */
		public double[][] d;
		
		/** 成本方程斜率。 */
		public double a;
		
		
		/** 成本方程截距。 */
		public double b;
		
		public ConstraintData(double[][] _d, double[] _q, double _a, double _b){
			d = _d;
			q = _q;
			a = _a;
			b = _b;
		}
	}
	

	public MinimizeFuelCost(ConstraintData _dat, Reference _ref, double _coef){
		super(_ref, false, _coef);
		dat =  _dat;
		fq = new double[ref.seq.length];
		bd = new double[ref.seq.length];
		name = "FUEL_COST_OBJECTIVE";
		update();
	}
	

	public MinimizeFuelCost(ConstraintData _dat, double _coef){
		super(false, _coef);
		dat = _dat;
		fq = null;
		bd = null;
		name = "DISTANCE_OBJECTIVE";
	}
	

	//************************************************************************************************
	//*************************feasibility check functions--------------------------------------------
	//************************************************************************************************
	public int insert_feasible(int p, int id){
		return 0;
	}
	

	public int remove_feasible(int p){
		return 0;
	}
	

	public int cross_feasible(int p1, Base cnt, int p2){
		return 0;
	}
	
	public int exchange_feasible(int p1, Base cnt, int p2){
		return 0;
	}
	

	public int replace_feasible(int p, int id){
		return 0;
	}
	
	
	//************************************************************************************************
	//*************************cost computation functions---------------------------------------------
	//************************************************************************************************
	public double insert_cost(int p, int id){
		int pre = ref.seq[p - 1];
		int next = ref.seq[p];
		double dist_delta = dat.d[pre][id] + dat.d[id][next] - dat.d[pre][next];
		double to_depot = dat.d[id][next] + bd[p];
		return coef * ((dat.a * fq[p - 1] + dat.b) * dist_delta + dat.a * dat.q[id] * to_depot);
	}
	
	public double remove_cost(int p){
		int pre = ref.seq[p - 1];
		int id = ref.seq[p];
		int next = ref.seq[p + 1];
		double dist_delta = dat.d[pre][next] - dat.d[pre][id] - dat.d[id][next];
		return coef * ((dat.a * fq[p - 1] + dat.b) * dist_delta - dat.a * dat.q[id] * bd[p]);
	}
	
	public double exchange_cost(int p1, Base obj, int p2){
		MinimizeFuelCost md = (MinimizeFuelCost) obj;
		int pre1 = ref.seq[p1 - 1];
		int id1 = ref.seq[p1];
		int next1 = ref.seq[p1 + 1];
		int pre2 = md.ref.seq[p2 - 1];
		int id2 = md.ref.seq[p2];
		int next2 = md.ref.seq[p2 + 1];
		double dist_delta1 = dat.d[pre1][id2] + dat.d[id2][next1] 
				- dat.d[pre1][id1] - dat.d[id1][next1];
		double dist_delta2 = md.dat.d[pre2][id1] + md.dat.d[id1][next2]
				- md.dat.d[pre2][id2] - md.dat.d[id2][next2];
		
		return coef * ((dat.a * fq[p1 - 1] + dat.b) * dist_delta1 + dat.a * dat.q[id2] * (dat.d[id2][next1] + bd[p1 + 1]) - dat.a * dat.q[id1] *  bd[p1])
				+ md.coef * ((md.dat.a * md.fq[p2 - 1] + md.dat.b) * dist_delta2 + md.dat.a * md.dat.q[id1] * (md.dat.d[id1][next2] + md.bd[p2 + 1]) - md.dat.a * md.dat.q[id2] *  md.bd[p2]);
				
	}
	
	public double cross_cost(int p1, Base obj, int p2){
		MinimizeFuelCost md = (MinimizeFuelCost) obj;
		int id1 = ref.seq[p1];
		int next1 = ref.seq[p1 + 1];
		int id2 = md.ref.seq[p2];
		int next2 = md.ref.seq[p2 + 1];
		double dist_delta1 = dat.d[id1][next2] + md.bd[p2 + 1] - bd[p1];
		double dist_delta2 = dat.d[id2][next1] + bd[p1 + 1] - md.bd[p2];
		return coef * ((dat.a * fq[p1] + dat.b) * dist_delta1) + md.coef * ((md.dat.a * md.fq[p2] + md.dat.b) * dist_delta2);
	}
	
	public double replace_cost(int p, int id){
		int pre = ref.seq[p - 1];
		int cur = ref.seq[p];
		int next = ref.seq[p + 1];
		double dist_delta = dat.d[pre][id] + dat.d[id][next] - dat.d[pre][cur] - dat.d[cur][next];
		return coef * ((dat.a * fq[p - 1] + dat.b) * dist_delta + dat.a * dat.q[id] * (dat.d[id][next] + bd[p + 1]) - dat.a * dat.q[cur] * bd[p]);
	}
	
	//************************************************************************************************
	//*************************query functions-------------------------------------------------------
	//************************************************************************************************
	public boolean is_constraint(){
		return false;
	}
	
	//************************************************************************************************
	//*************************functions to change the status of the objective------------------------
	//************************************************************************************************
	public void update(){
		c = 0;
		fq[0] = dat.q[ref.seq[0]];
		for(int i = 0; i < ref.len - 1; i++){
			int id = ref.seq[i];
			int next = ref.seq[i + 1];
			if(i > 0)
				fq[i] = fq[i - 1] + dat.q[id];
			c += coef * (dat.a * fq[i] + dat.b) * dat.d[id][next];
		}
		bd[ref.len - 1] = 0;
		for(int i = ref.len - 2; i >= 0; i--){
			bd[i] = bd[i + 1] + dat.d[ref.seq[i]][ref.seq[i + 1]];
		}
	}

	public Constraint copy(Reference _ref){
		return new MinimizeFuelCost(new ConstraintData(dat.d, dat.q, dat.a, dat.b), _ref, coef);
	}
	
	
}
