import gzip
import tarfile
import os

import gurobipy as gp
from gurobipy import GRB
 
class Data():
    ''' Read data. '''
    def __init__(self): 
        # Description of tar file containing instances.
        self.name_tar_file = "BOBILib/general-bilevel.tar.gz"
        self.folders_to_read = ["general-bilevel/pure-integer/denegre/",\
                                "general-bilevel/pure-integer/zhang/",\
                                "general-bilevel/mixed-integer/xuwang/"] 
        self.init_name_instance = {"general-bilevel/pure-integer/denegre/": "denegre",\
                                   "general-bilevel/pure-integer/zhang/": "zhang",\
                                   "general-bilevel/mixed-integer/xuwang/": "xuwang" } 

        # Reading instances at defined folders.
        with tarfile.open(self.name_tar_file, "r:gz") as tar:
            for member in tar.getmembers():
                # Skip directories.
                if member.isdir():
                    continue
                
                # Check folder and extension.
                if any(member.name.startswith(folder) for folder in self.folders_to_read) and member.name.endswith(".mps.gz"):
                    
                    # Base name of current instance.
                    name_instance = member.name.rsplit(".mps.gz", 1)[0]  

                    # Write temporary files for current instance.
                    self.writeTempFiles(tar,member,name_instance)

                    # Read temporary files for current instance.
                    vars_info, constrs_info, follower_info = self.readTempFiles()
                    
                    # Define instance information.
                    leader_vars, follower_vars, leader_constrs, follower_constrs, leader_lb, leader_ub, follower_lb, follower_ub  = self.defineInstance(vars_info,constrs_info,follower_info)

                    # Write instance file.
                    folder_match = next((folder for folder in self.folders_to_read if member.name.startswith(folder)), None)
                    name_instance = self.init_name_instance[folder_match] + "-" + os.path.basename(name_instance)
                    self.writeInstanceFile(name_instance,leader_vars,follower_vars,leader_constrs,follower_constrs,leader_lb,leader_ub,follower_lb,follower_ub)

                    # Remove temporary files.
                    os.remove("temp.aux")
                    os.remove("temp.mps")

    ''' Method for reading current instance .mps and .aux files and copying them as temporary files. '''
    def writeTempFiles(self,tar,member,name_instance):
        # Extract and copy .mps file.
        print(f"MPS file: {member.name}")
        f = tar.extractfile(member)
        with gzip.open(f, "rt") as mps_file:
            with open("temp.mps", "w") as temp:
                temp.write(mps_file.read())
        
        # Try to find the corresponding .aux file in the tar
        aux_member_name = name_instance + ".aux"
        try:
            aux_member = tar.getmember(aux_member_name)
            aux_f = tar.extractfile(aux_member)
            if aux_f is not None:
                print(f"AUX file: {aux_member.name}")
                with open("temp.aux", "wb") as temp_aux:
                    temp_aux.write(aux_f.read())
        except KeyError:
            raise RuntimeError(f"Failed to find file: {aux_member_name}")

    ''' Method for reading information on .mps and .aux temporary files. '''
    def readTempFiles(self):
        vars_info, constrs_info = self.readTempMpsFile()
        follower_info = self.readTempAuxFile()
        return vars_info, constrs_info, follower_info
    
    ''' Method for reading information on .mps temporary file. '''
    def readTempMpsFile(self):
        model = gp.read("temp.mps") # using gurobi to read file.
        vars_info = {v.VarName : {"obj_leader": v.Obj, "obj_follower": 0.0, "lb": v.LB, "ub": v.UB,"type": v.VType} for v in model.getVars()}
        constrs_info = {c.ConstrName : {"sense": c.Sense, "rhs": c.RHS, "coeffs": {v.VarName: model.getCoeff(c, v) for v in model.getVars() if model.getCoeff(c, v) != 0}} for c in model.getConstrs()}
        obj_info = model.ModelSense

        # Always define a minimization problem
        if obj_info == GRB.MAXIMIZE:
            for vname, vinfo in vars_info.items():
                vinfo["obj_leader"] = -vinfo["obj_leader"]

        return vars_info, constrs_info

    ''' Method for reading information on .aux temporary file. '''
    def readTempAuxFile(self):
        # Data .aux file.
        data = {"num_vars": 0, "num_constrs": 0, "vars": {}, "constrs": [], "name": "", "mps_file": ""}

        # Auxiliary function to read .aux file.
        def parse_aux_line(line, section):
            if section == "num_vars":
                data["num_vars"] = int(line)
            elif section == "num_constrs":
                data["num_constrs"] = int(line)
            elif section == "vars":
                parts = line.split()
                var_name = parts[0]
                value = float(parts[1])
                data["vars"][var_name] = value
            elif section == "constrs":
                data["constrs"].append(line)
            elif section == "name":
                data["name"] = line
            elif section == "mps_file":
                data["mps_file"] = line

        # Read .aux file.
        with open("temp.aux", "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        section = None
        for line in lines:
            if line.startswith("@NUMVARS"):
                section = "num_vars"
            elif line.startswith("@NUMCONSTRS"):
                section = "num_constrs"
            elif line.startswith("@VARSBEGIN"):
                section = "vars"
            elif line.startswith("@VARSEND"):
                section = None
            elif line.startswith("@CONSTRSBEGIN"):
                section = "constrs"
            elif line.startswith("@CONSTRSEND"):
                section = None
            elif line.startswith("@NAME"):
                section = "name"
            elif line.startswith("@MPS"):
                section = "mps_file"
            elif section is not None:
                parse_aux_line(line, section)
        
        return data

    ''' Method defining current instance information. '''
    def defineInstance(self,vars,constrs,data_follower):
        
        ''' Auxiliar function returning value considered for follower variable to modify coupling constraints.'''
        def follower_guess(var):
            lb = var["lb"]
            ub = var["ub"]
            if lb <= -GRB.INFINITY and ub >= GRB.INFINITY:
                guess = 0.0
            elif lb <= -GRB.INFINITY:
                guess = ub
            elif ub >= GRB.INFINITY:
                guess = lb
            else:
                guess = (lb + ub) / 2.0
            return guess

        # Create dictionaries with leader and follower variables.
        leader_vars = {}
        follower_vars = {}
        for var_name, info in vars.items():
            if var_name in data_follower["vars"]: # This is a follower variable.
                # Transform variable into continuous to obtain linear follower problem.
                info["type"] = 'C'  
                # Define follower objective.                               
                info["obj_follower"] = data_follower["vars"][var_name] 
                follower_vars[var_name] = info
            else: # This is a leader variable.
                leader_vars[var_name] = info

        # Create dictionaries with leader and follower constraints.
        leader_constrs = {}
        follower_constrs = {}
        for constr_name, info in constrs.items():
            if constr_name in data_follower["constrs"]: # This is a follower constraint.
                follower_constrs[constr_name] = info
            else: # This is a leader constraint.
                # Remove coupling terms
                if any(var in follower_vars for var in info["coeffs"]):
                    for var in list(info["coeffs"]):
                        if var in follower_vars:
                            y_guess = follower_guess(follower_vars[var])
                            info["rhs"] -= info["coeffs"][var] * y_guess
                            del info["coeffs"][var]
                leader_constrs[constr_name] = info

        leader_lb = self.getLeaderLowerBound(leader_vars,follower_vars,leader_constrs,follower_constrs)
        leader_ub = self.getLeaderUpperBound(leader_vars,follower_vars,leader_constrs,follower_constrs)
        
        follower_lb = self.getFollowerLowerBound(leader_vars,follower_vars,leader_constrs,follower_constrs)
        follower_ub = self.getFollowerUpperBound(leader_vars,follower_vars,leader_constrs,follower_constrs)
        
        return leader_vars, follower_vars, leader_constrs, follower_constrs, leader_lb, leader_ub, follower_lb, follower_ub 

    ''' Compute lower bound on follower problem from current instance. Assuming follower is a minimization problem. '''
    def getFollowerLowerBound(self,leader_vars,follower_vars,leader_constrs,follower_constrs):
        
        # Solve gurobi model to obtain lower bound on follower objective.
        # We want the leader decision that yields the mininum optimal value for the follower objective.
        flb_model = gp.Model()

        x = flb_model.addVars([name for name, _ in leader_vars.items()],
                    lb=[info["lb"] for _, info in leader_vars.items()],
                    ub=[info["ub"] for _, info in leader_vars.items()],
                    obj=[info["obj_follower"] for _, info in leader_vars.items()],
                    vtype=[info["type"] for _, info in leader_vars.items()],
                    name=[name for name, _ in leader_vars.items()])

        y = flb_model.addVars([name for name, _ in follower_vars.items()],
                    lb=[info["lb"] for _, info in follower_vars.items()],
                    ub=[info["ub"] for _, info in follower_vars.items()],
                    obj=[info["obj_follower"] for _, info in follower_vars.items()],
                    vtype=[info["type"] for _, info in follower_vars.items()],
                    name=[name for name, _ in follower_vars.items()])
        
        for cname, cinfo in leader_constrs.items():
            if cinfo["sense"] == "<":
                flb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) <= cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                flb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) >= cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                flb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) == cinfo["rhs"], cname)

        for cname, cinfo in follower_constrs.items():
            if cinfo["sense"] == "<":
                flb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) <= cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                flb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) >= cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                flb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) == cinfo["rhs"], cname)

        # Defining objective.
        flb_model.modelSense = GRB.MINIMIZE

        # Solve follower lower bound problem.
        flb_model.optimize()

        # Obtain lower bound
        follower_lb = -GRB.INFINITY
        if flb_model.Status == GRB.OPTIMAL:
            follower_lb = flb_model.ObjVal
        
        return follower_lb

    ''' Compute upper bound on follower problem from current instance. Assuming follower is a minimization problem. '''
    def getFollowerUpperBound(self,leader_vars,follower_vars,leader_constrs,follower_constrs):
        
        # Solve gurobi model to obtain upper bound on follower objective.
        # We want the leader decision that yields the maximum optimal value for the follower objective. (max-min problem)
        fub_model = gp.Model()

        # Leader variables and constraints.
        x = fub_model.addVars([name for name, _ in leader_vars.items()],
                    lb=[info["lb"] for _, info in leader_vars.items()],
                    ub=[info["ub"] for _, info in leader_vars.items()],
                    obj=[info["obj_follower"] for _, info in leader_vars.items()],
                    vtype=[info["type"] for _, info in leader_vars.items()],
                    name=[name for name, _ in leader_vars.items()])
        
        for cname, cinfo in leader_constrs.items():
            if cinfo["sense"] == "<":
                fub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) <= cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                fub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) >= cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                fub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) == cinfo["rhs"], cname)

        # Follower variables and constraints.
        y = fub_model.addVars([name for name, _ in follower_vars.items()],
                    lb=[info["lb"] for _, info in follower_vars.items()],
                    ub=[info["ub"] for _, info in follower_vars.items()],
                    obj=[info["obj_follower"] for _, info in follower_vars.items()],
                    vtype=[info["type"] for _, info in follower_vars.items()],
                    name=[name for name, _ in follower_vars.items()])

        slack = {}
        for cname, cinfo in follower_constrs.items():
            if cinfo["sense"] == "<":
                slack[cname] = fub_model.addVar(lb=0, name=f"slack_{cname}")
                fub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) + slack[cname] == cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                slack[cname] = fub_model.addVar(lb=0, name=f"slack_{cname}")
                fub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) - slack[cname] == cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                fub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) == cinfo["rhs"], cname)

        slack_lb = {}
        slack_ub = {}
        for vname, vinfo in follower_vars.items():
            if vinfo["lb"] > -GRB.INFINITY:
                slack_lb[vname] = fub_model.addVar(lb=0, name=f"slack_lb_{vname}")
                fub_model.addConstr(y[vname] - slack_lb[vname] == vinfo["lb"])
            if vinfo["ub"] < GRB.INFINITY:
                slack_ub[vname] = fub_model.addVar(lb=0, name=f"slack_ub_{vname}")
                fub_model.addConstr(y[vname] + slack_ub[vname] == vinfo["ub"])


        # Dual variables for follower constraints.
        u = {}
        for cname, cinfo in follower_constrs.items():
            if cinfo["sense"] == "<":
                u[cname] = fub_model.addVar(lb=-GRB.INFINITY, ub=0, name=f"u_{cname}")
            elif cinfo["sense"] == ">":
                u[cname] = fub_model.addVar(lb=0, ub=GRB.INFINITY, name=f"u_{cname}")
            elif cinfo["sense"] == "=":
                u[cname] = fub_model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"u_{cname}")
        
        # Dual variables from lower and upper bounds on follower variables.
        u_lb = {}
        u_ub = {}
        for vname, vinfo in follower_vars.items():
            if vinfo["lb"] > -GRB.INFINITY:
                u_lb[vname] = fub_model.addVar(lb=0, ub=GRB.INFINITY, name=f"u_lb_{vname}")
            if vinfo["ub"] < GRB.INFINITY:
                u_ub[vname] = fub_model.addVar(lb=0, ub=GRB.INFINITY, name=f"u_ub_{vname}")

        # Dual constraints for follower variables
        for vname, vinfo in follower_vars.items():
            expr = gp.LinExpr()
            for cname, cinfo in follower_constrs.items():
                coef = cinfo["coeffs"].get(vname, 0.0)
                expr += coef * u[cname]
            
            if vinfo["lb"] <= -GRB.INFINITY and vinfo["ub"] >= GRB.INFINITY:
                fub_model.addConstr(expr == vinfo["obj_follower"], name=f"dual_constr_{vname}")
            elif vinfo["lb"] <= -GRB.INFINITY:
                fub_model.addConstr(expr - u_ub[vname] == vinfo["obj_follower"], name=f"dual_constr_{vname}")
            elif vinfo["ub"] >= GRB.INFINITY:
                fub_model.addConstr(expr + u_lb[vname] == vinfo["obj_follower"], name=f"dual_constr_{vname}")
            else:
                fub_model.addConstr(expr + u_lb[vname] - u_ub[vname] == vinfo["obj_follower"], name=f"dual_constr_{vname}")
        
        # Complemenary slackness.
        for cname, cinfo in follower_constrs.items():
            if cinfo["sense"] != "=":
                fub_model.addSOS(GRB.SOS_TYPE1, [u[cname], slack[cname]])
        
        for vname, vinfo in follower_vars.items():
            if vinfo["lb"] > -GRB.INFINITY:
                fub_model.addSOS(GRB.SOS_TYPE1, [u_lb[vname], slack_lb[vname]])
            if vinfo["ub"] < GRB.INFINITY:
                fub_model.addSOS(GRB.SOS_TYPE1, [u_ub[vname], slack_ub[vname]])
        
        # Defining objective
        fub_model.modelSense = GRB.MAXIMIZE

        # Solving follower upper bound problem.
        fub_model.optimize()

        # Obtain upper bound
        follower_ub = GRB.INFINITY
        if fub_model.Status == GRB.OPTIMAL:
            follower_ub = fub_model.ObjVal
        
        return follower_ub

    ''' Compute any lower bound on leader problem (considering follower variable only on the objective) from current instance. Assuming is a minimization problem. '''
    def getLeaderLowerBound(self,leader_vars,follower_vars,leader_constrs,follower_constrs):
        # Solve gurobi model to obtain lower bound on leader objective.
        # Relaxing optimality constraint on follower problem and minimizing it.
        llb_model = gp.Model()

        # Time limit of 1 hour
        llb_model.setParam('timelimit',3600)

        x = llb_model.addVars([name for name, _ in leader_vars.items()],
                    lb=[info["lb"] for _, info in leader_vars.items()],
                    ub=[info["ub"] for _, info in leader_vars.items()],
                    # obj=[info["obj_leader"] for _, info in leader_vars.items()],
                    obj=[0.0 for _, info in leader_vars.items()],
                    vtype=[info["type"] for _, info in leader_vars.items()],
                    name=[name for name, _ in leader_vars.items()])

        y = llb_model.addVars([name for name, _ in follower_vars.items()],
                    lb=[info["lb"] for _, info in follower_vars.items()],
                    ub=[info["ub"] for _, info in follower_vars.items()],
                    obj=[info["obj_leader"] for _, info in follower_vars.items()],
                    vtype=[info["type"] for _, info in follower_vars.items()],
                    name=[name for name, _ in follower_vars.items()])
        
        for cname, cinfo in leader_constrs.items():
            if cinfo["sense"] == "<":
                llb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) <= cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                llb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) >= cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                llb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) == cinfo["rhs"], cname)

        for cname, cinfo in follower_constrs.items():
            if cinfo["sense"] == "<":
                llb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) <= cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                llb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) >= cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                llb_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) == cinfo["rhs"], cname)

        # Defining objective.
        llb_model.modelSense = GRB.MINIMIZE

        # Solve follower lower bound problem.
        llb_model.optimize()

        # Obtain lower bound
        leader_lb = -GRB.INFINITY
        if llb_model.Status == GRB.OPTIMAL or llb_model.Status == GRB.TIME_LIMIT:
            leader_lb = llb_model.ObjVal
        
        return leader_lb 

    ''' Compute any upper bound on leader problem (considering follower variable only on the objective) from current instance. Assuming is a minimization problem. '''
    def getLeaderUpperBound(self,leader_vars,follower_vars,leader_constrs,follower_constrs):
        # Solve gurobi model to obtain upper bound on leader objective.
        # Relaxing optimality constraint on follower problem and maximixing leader objective.
        lub_model = gp.Model()

        # Time limit of 1 hour
        lub_model.setParam('timelimit',3600)

        x = lub_model.addVars([name for name, _ in leader_vars.items()],
                    lb=[info["lb"] for _, info in leader_vars.items()],
                    ub=[info["ub"] for _, info in leader_vars.items()],
                    # obj=[info["obj_leader"] for _, info in leader_vars.items()],
                    obj=[0.0 for _, info in leader_vars.items()],
                    vtype=[info["type"] for _, info in leader_vars.items()],
                    name=[name for name, _ in leader_vars.items()])

        y = lub_model.addVars([name for name, _ in follower_vars.items()],
                    lb=[info["lb"] for _, info in follower_vars.items()],
                    ub=[info["ub"] for _, info in follower_vars.items()],
                    obj=[info["obj_leader"] for _, info in follower_vars.items()],
                    vtype=[info["type"] for _, info in follower_vars.items()],
                    name=[name for name, _ in follower_vars.items()])
        
        for cname, cinfo in leader_constrs.items():
            if cinfo["sense"] == "<":
                lub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) <= cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                lub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) >= cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                lub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) == cinfo["rhs"], cname)

        for cname, cinfo in follower_constrs.items():
            if cinfo["sense"] == "<":
                lub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) <= cinfo["rhs"], cname)
            elif cinfo["sense"] == ">":
                lub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) >= cinfo["rhs"], cname)
            elif cinfo["sense"] == "=":
                lub_model.addConstr(gp.quicksum(cinfo["coeffs"].get(vname,0.0)*x[vname] for vname in leader_vars.keys()) + gp.quicksum(cinfo["coeffs"].get(vname,0.0)*y[vname] for vname in follower_vars.keys()) == cinfo["rhs"], cname)

        # Defining objective.
        lub_model.modelSense = GRB.MAXIMIZE

        # Solve follower lower bound problem.
        lub_model.optimize()

        # Obtain lower bound
        leader_ub = GRB.INFINITY
        if lub_model.Status == GRB.OPTIMAL or lub_model.Status == GRB.TIME_LIMIT:
            leader_ub = lub_model.ObjVal
        
        return leader_ub 

    ''' Write current instance file. '''
    def writeInstanceFile(self,name_instance,leader_vars,follower_vars,leader_constrs,follower_constrs,leader_lb,leader_ub,follower_lb,follower_ub):
        name = 'instances/' + name_instance + '.txt'
        with open(name, 'w') as f:
            f.write( "%d %d %.3f %.3f\n" % (len(leader_vars),len(leader_constrs),leader_lb,leader_ub))
            f.write( "%d %d %.3f %.3f\n" % (len(follower_vars),len(follower_constrs),follower_lb,follower_ub))
            for vname, vinfo in leader_vars.items():
                f.write( "%s %.3f %.3f %.3f %s\n" % (vname,vinfo["lb"],vinfo["ub"],vinfo["obj_leader"],vinfo["type"]))
            for vname, vinfo in follower_vars.items():
                f.write( "%s %.3f %.3f %.3f %.3f %s\n" % (vname,vinfo["lb"],vinfo["ub"],vinfo["obj_leader"],vinfo["obj_follower"],vinfo["type"]))
            for cname, cinfo in leader_constrs.items():
                f.write("%s %s %.3f " % (cname,cinfo["sense"],cinfo["rhs"]))
                for vname, vinfo in leader_vars.items():
                    f.write("%.3f " % (cinfo["coeffs"].get(vname,0.0)))
                f.write("\n")
            for cname, cinfo in follower_constrs.items():
                f.write("%s %s %.3f " % (cname,cinfo["sense"],cinfo["rhs"]))
                for vname, vinfo in leader_vars.items():
                    f.write("%.3f " % (cinfo["coeffs"].get(vname,0.0)))
                for vname, vinfo in follower_vars.items():
                    f.write("%.3f " % (cinfo["coeffs"].get(vname,0.0)))
                f.write("\n")

def main():
    data = Data()
    
if __name__ == "__main__":
    main()