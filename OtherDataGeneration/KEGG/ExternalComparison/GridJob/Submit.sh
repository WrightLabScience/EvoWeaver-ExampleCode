universe = vanilla
executable = Wrapper.sh
log = LogFiles/out.$(Process).log
output = LogFiles/out.$(Process).out
error = LogFiles/out.$(Process).err

# Kernel version comes from major.minor.patch in the format:
# major * 10000 + minor * 1000 + patch
# SynExtend comes from r-base which uses debian:latest, which
# as of 20190722 is "buster" from 4.19.105
# as per advice, 31000 will keep jobs away from RHEL 6

requirements = Arch == "X86_64" && HAS_SINGULARITY == True && OSG_HOST_KERNEL_VERSION >= 31000
request_cpus = 1
request_memory = 2GB
request_disk = 8GB

+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/npcooley/synextend:slim.1.20.0"

# IF using a docker container, only send the RScript and any associated data files
transfer_input_files = JobScript.R, \
                        # inputfiles/$(FileName), \
                        # stash:///osgconnct/public/ahl/$(FileName), \
                        # osdf:///ospool/ap20/data/ahl/$(FileName)


# Send jobs to Held state on failure up to 4 times, a 5th attempt is allowed but afterwards it should
# register as completed
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0) && \
               (NumJobStarts < 5)

# release a job that has been held, after it has been held for 10 minutes, up to a maximum of 5 retries
#periodic_release = (NumJobStarts < 10) && \
#                   (((HoldReasonCode == 12) && (HoldReasonSubCode == 0 )) || (HoldReasonCode == 3)) && \
#                   ((CurrentTime - EnteredCurrentStatus) > 600)

# don't overwhelm the submit node
# max_idle = 1000
max_materialize = 10000

# arguments = $(Process) $(FileName)
arguments = $(Process) $(Organism)

queue Process, Organism from JobMap.txt
#queue 5
