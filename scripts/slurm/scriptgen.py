import re

pipe_pattern = re.compile("[|]")
def escape_pipes_filenames(path):
    return re.sub(pipe_pattern, "\|", path)

class SlurmScriptGenerator(object):
    def __init__(
            self,
            jobname="somejerb",
            nodes=1,
            tasks_per_node=1,
            cpus_per_task=1,
            mem_per_cpu=2,
            time=10,
            account="tyjames1",
            partition="standard"
            ):

        self.jobname = jobname
        self.header = {
                "#SBATCH --job-name": jobname,
                "#SBATCH --mail-type": "BEGIN,END",
                "#SBATCH --nodes": nodes,
                "#SBATCH --ntasks-per-node": tasks_per_node,
                "#SBATCH --cpus-per-task": cpus_per_task,
                "#SBATCH --mem-per-cpu": f"{mem_per_cpu}g",
                "#SBATCH --time": f"{time}:00:00",
                "#SBATCH --account": account,
                "#SBATCH --partition": partition
                }

        self.commands = list()

    def add_command(self, cmd) -> None:
        if isinstance(cmd, list):
            self.commands.append(" ".join(cmd).strip())
        elif isinstance(cmd, str):
            self.commands.append(cmd.strip())
        else:
            raise TypeError("Command must be list or str.")

        return None

    def __repr__(self):
        nl = "\n"
        stringify = [f"{k}={v}" for k,v in self.header.items()]
        str_header = f"{nl.join(stringify)}"
        return f"#!/bin/bash{nl}{str_header}{nl*2}{nl.join(self.commands)}"

    def write(self):
        with open(f"{self.jobname}.sh", 'w') as f:
            f.write(self.__repr__())
