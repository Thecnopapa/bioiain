import os, sys, json, asyncio, time, threading, ctypes, psutil
from ..utilities.logging import log


#print("START")
log("header", "Importing parallel utils...")
cpu_count = psutil.cpu_count(logical=False)
logical_cpu = psutil.cpu_count(logical=True)
gpu_count = 0
is_cluster=False
use_max = False
use_logical = False
if use_logical:
    cpu_count = logical_cpu
if os.environ.get("SLURM_JOB_ID", None) is not None:
    log(1, "SLURM manager detected, maximizing CPUS...")
    #for k, v  in os.environ.items():
    #    if k.startswith("SLURM"):
    #        print(k, v)
    cpu_count = int(os.environ["SLURM_CPUS_ON_NODE"])
    try:
        gpu_count = int(os.environ["SLURM_GPUS_ON_NODE"])
    except: pass
    use_max = True
    is_cluster = True


if use_max or cpu_count <= 1:
    avail_cpus = cpu_count
else:
    avail_cpus = cpu_count -1
log(1, f"Available CPUs: {avail_cpus}/{cpu_count}, using max: {use_max}, using_logical: {use_logical}")


pools = []

def mem_usage(as_dict=False):
    if as_dict:
        return psutil.virtual_memory()._asdict()
    else:
        return psutil.virtual_memory().available * 100 / psutil.virtual_memory().total




def _indefinite_mem_log():
    print("MEMORY LOGGING STARTED")

    open("./mem_log.txt", "w")
    while True:
        with open("./mem_log.txt", "a") as f:
            d = mem_usage(as_dict=True)
            f.write(f"{d['percent']:2.2f}% ({d['active']/1000000000:2.2f} GB)\n")
        time.sleep(1)



def mem_log():
    pool = ThreadPool()
    pool.add(_indefinite_mem_log)
    print(pool)
    pools.append(pool)
    pool.start(wait=False)





def split_iterable(iterable, n_parts:int|str="auto"):

    if n_parts == "auto":
        if use_max:
            n_parts = cpu_count
        else:
            n_parts = cpu_count - 1
    elif n_parts == "max":
        n_parts = cpu_count
    elif n_parts == "double":
        n_parts = cpu_count*2

    assert type(n_parts) == int
    if n_parts <= 1:
        return [iterable]

    l = len(iterable)
    log("header", "Splitting iterable of length:", l)
    log(1, "Number of parts:", n_parts,)
    part_size = l//n_parts
    last_part_size = l%n_parts
    log(1, "Size of parts:", part_size, f"({last_part_size})" )
    print(part_size, last_part_size)
    out = []
    t = type(iterable)

    for n in range(n_parts):
        start = part_size*n
        if n == n_parts -1:
            if last_part_size == 0:
                continue
            end = start + last_part_size
        else:
            end = start + part_size
        if t in (list, tuple, str):
            out.append([e for e in iterable[start:end]])
        elif t == dict:
            out.append({k:v for k,v in iterable.items()[start:end]})
        else:
            log("error", "Unrecognised iterable type:", t)
            return iterable
    return out






class AsyncPool(object):
    """
    To be used exclusively with async functions
    """
    def __init__(self, n_workers="auto"):
        self.n_workers = n_workers
        self.tasks = {}
        self.current_pool = None
        self.current_keys = None
        self.current_tasks = None
        self.running = False
        self.start_time = None
        self.end_time = None

    def __len__(self):
        return len(self.tasks)

    def __repr__(self):
        pending = 0
        running = 0
        done = 0
        errors = 0
        for task in self.tasks.values():
            if task["status"] == "running":
                running += 1
            elif task["status"] == "done":
                done += 1
            elif task["status"] == "error":
                errors += 1
            elif task["status"] == "pending":
                pending += 1
        return f"<AsyncPool: pending:{pending} running: {running} done:{done} errors:{errors}>"

    def add(self, awaitable, task_id=None):
        if task_id is None:
            task_id = len(self.tasks.keys())
        task_id = str(task_id)
        if task_id in self.tasks.keys():
            raise Exception(f"Task {task_id} already exists")
        self.tasks[task_id] = {"awaitable": awaitable, "status": "pending"}
        return self.tasks[task_id]

    def __add__(self, awaitable, task_id=None):
        return self.add(awaitable, task_id=task_id)



    def info(self):
        print(repr(self)[:-1], end="\n")
        for k, task in self.tasks.items():
            print(f" - Task {k} ({task['status']}): {task['awaitable'].__name__}", end="")
            if "return" in task.keys():
                print(f" --> ({type(task['return']).__name__})", end="")
                if type(task["return"]) in [list, dict, tuple]:
                    print(" of length:", len(task["return"]), end="")
                else:
                    print(":",task["return"], end="")
            print("")
        print(">")


    async def _run(self, wait=False, **kwargs):
        if self.running:
            log("warning", "AsyncPool.run(): Already running")
            raise Exception(f"Already running")
        self.running = True
        self.current_keys = []
        self.current_tasks = []
        for k, v in self.tasks.items():
            if v["status"] != "pending":
                continue
            self.current_keys.append(k)
            self.current_tasks.append(v["awaitable"])
            self.tasks[k]["status"] = "running"
        n_tasks = len(self.current_tasks)
        log("header", f"AsyncPool: Running {n_tasks} tasks (wait={wait})")
        self.current_pool = asyncio.gather(*self.current_tasks, return_exceptions=True)
        if True:
            return await self._await(**kwargs)
        else:
            return self._await(**kwargs)




    async def _await(self, raise_errors=True, return_dict=False, **kwargs):
        if self.current_pool is None:
            log("warning", "AsyncPool._await(): No running pool")
            return None
        ret = await self.current_pool
        errors = 0
        ok = 0
        for k, rv in zip(self.current_keys, ret):
            self.tasks[k]["return"] = rv
            if isinstance(rv, Exception):
                self.tasks[k]["status"] = "error"
                errors += 1
                log("error", f"in task: {k}: {rv}")
                if raise_errors:
                    raise rv
            else:
                self.tasks[k]["status"] = "done"
                ok += 1
        log("header", f"AsyncPool: Finished {ok+errors} tasks ({errors} errors)")
        self.running = False
        self.current_keys = None
        self.current_tasks = None
        self.current_pool = None

        if return_dict:
            return self.tasks
        return ret




    def start(self, wait=False, **kwargs):
        asyncio.run(self._run(wait=wait, **kwargs))
        return self

    def without_errors(self, raise_errors=False, **kwargs):
        return self.start(raise_errors=raise_errors, **kwargs)

    def _deprecated_wait(self, wait=True, **kwargs):
        return self.start(wait=wait, **kwargs)

    def wait(self, **kwargs):
        asyncio.run(self._await(**kwargs))
        return self


    def get_return(self, key):
        return self.tasks[key]["return"]





class ThreadPool(object):
    def __init__(self):
        self.threads = {}
        self.running = False
        self.context = None
        self.returns = {}
        pools.append(self)



    class Thread(threading.Thread):
        def __init__(self, *args, ret=None, **kwargs):
            super().__init__(*args, **kwargs)
            self.ret = ret

            self.error = False


        class ThreadKilled(BaseException):
            def __init__(self, *args, **kwargs):
                super().__init__(self, *args, **kwargs)
                log("warning", f"Thread Killed:")



        def run(self, *args, **kwargs):
            try:
                r = super().run(*args, **kwargs)
                self.ret = r
                return r
            except:
                print("ERRRRORRRR")
                self.error = True
                raise
                return None

        def get_id(self):

            if hasattr(self, '_thread_id'):
                print("TREAD_ID:", self._thread_id)
                return self._thread_id
            for id, thread in threading._active.items():
                if thread is self:
                    return id

        def kill(self):
            ptssae = ctypes.pythonapi.PyThreadState_SetAsyncExc
            ptssae.argtypes = (ctypes.c_ulong, ctypes.py_object)
            ptssae.restype = ctypes.c_int
            thread_id = self.get_id()
            log("warning", "Trying to kill thread:", thread_id)
            res = ptssae(thread_id, ctypes.py_object(self.ThreadKilled))

            return res


    def add(self, fun, *args, **kwargs):
        t = self.Thread(target=fun, *args, context=self.context, **kwargs)
        self.threads[len(self.threads)] = {"thread": t, "fun": fun, "status": "pending", "ret":t.ret, "name":None}
        return t

    def start(self, wait=False, **kwargs):
        pending_threads = {k: v for k, v in self.threads.items() if v["status"] == "pending"}
        for k, t in pending_threads.items():
            t["thread"].name = "Thread {}".format(k)
            t["name"] = t["thread"].name
            t["thread"].start()
            t["status"] = "running"
        log("header", f"ThreadPool: Running {len(pending_threads)} tasks (wait={wait})")
        if wait:
            self._await(**kwargs)


    def _await(self, **kwargs):
        running_threads = {k:v for k,v in self.threads.items() if v["status"] == "running"}
        ok = 0
        errors = 0
        for k, t in running_threads.items():
            t["thread"].join()

            self.returns[t["name"]] = t["ret"]
            if t["thread"].error:
                t["status"] = "error"
                errors += 1
            else:
                t["status"] = "done"
                ok += 1
        log("header", f"Threadpool: Finished {ok + errors} tasks ({errors} errors), returning:")
        [log(1, f"{k}:", r) for k, r in self.returns.items()]
        pools.remove(self)
        return self.returns.values()

    def terminate(self):
        running_threads = {k:v for k,v in self.threads.items() if v["status"] == "running"}
        for k, t in running_threads.items():
            #print(t["thread"])
            t["thread"].kill()
        self._await()



def end_pools():
    if len(pools) > 0:
        log("start", "ENDING POOLS")
        #print(pools)
        for pool in pools:
            #print(pool)
            pool.terminate()

        log("end", "POOLS TERMINATED")


