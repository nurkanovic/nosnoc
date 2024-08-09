from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import json
# TODO portableize this
nosnoc_root = Path("/home/anton/syscop/software/nosnoc")

def gen_cmake(env, solver_data):
    template = env.get_template("CMakeLists.txt")
    content = template.render(solver_data)

    with open("CMakeLists.txt", mode="w", encoding="utf-8") as cmakelists:
        cmakelists.write(content)
        #print(f"{content}")

def gen_HS_cpp(env, solver_data):
    template = env.get_template("HomotopySolver.cpp")
    content = template.render(solver_data)

    with open("HomotopySolver.cpp", mode="w", encoding="utf-8") as cmakelists:
        cmakelists.write(content)
        #print(f"{content}")

def gen_HS_hpp(env, solver_data):
    template = env.get_template("HomotopySolver.hpp")
    content = template.render(solver_data)

    with open("HomotopySolver.hpp", mode="w", encoding="utf-8") as cmakelists:
        cmakelists.write(content)
        #print(f"{content}")


if __name__=="__main__":
    print(str(Path(nosnoc_root,'codegen_templates')))
    env = Environment(loader=FileSystemLoader(str(Path(nosnoc_root,'codegen_templates'))))
    
    with open("solver.json") as f:
        solver_data = json.load(f)

    if not isinstance(solver_data["nlp_lbw"], list):
        solver_data["nlp_lbw"] = [solver_data["nlp_lbw"]];
    if not isinstance(solver_data["nlp_ubw"], list):
        solver_data["nlp_ubw"] = [solver_data["nlp_ubw"]];
    if not isinstance(solver_data["nlp_lbg"], list):
        solver_data["nlp_lbg"] = [solver_data["nlp_lbg"]];
    if not isinstance(solver_data["nlp_ubg"], list):
        solver_data["nlp_ubg"] = [solver_data["nlp_ubg"]];
    if not isinstance(solver_data["nlp_p0"], list):
        solver_data["nlp_p0"] = [solver_data["nlp_p0"]];
    if not isinstance(solver_data["nlp_x0"], list):
        solver_data["nlp_x0"] = [solver_data["nlp_x0"]];

    gen_cmake(env, solver_data)
    gen_HS_cpp(env, solver_data)
    gen_HS_hpp(env, solver_data)
    
    
