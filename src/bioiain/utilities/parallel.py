import os, sys, json, asyncio, time, threading, ctypes, psutil, platform
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
log(2,f"Machine: {platform.processor()} ({platform.machine()})")


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





def split_iterable(iterable, n_parts:int|str="auto") -> list:

    if n_parts == "auto":
        n_parts = avail_cpus
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
    if l <= n_parts:
        n = l
        part_size = 1
        last_part_size = 0
    else:
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





class ThreadPool(object):
    def __init__(self):
        self.threads = {}
        self.running = False
        self.context = None
        self.returns = {}
        pools.append(self)



    class Thread(threading.Thread):
        def __init__(self, *args, target=None, ret=None, **kwargs):
            super().__init__(target=target, args=args, kwargs=kwargs)
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
        t = self.Thread(*args, target=fun, **kwargs)
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
        for k, t in running_threads.items():
            log(1, f"{k} ({t['status']}): {t['ret']}")

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


