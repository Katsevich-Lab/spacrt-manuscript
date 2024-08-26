# Interval in seconds for checking the job status
CHECK_INTERVAL=600

# Job name patterns to search for
searchPatterns="GCM|spaCRT|score|power|sceptre"

sleep $CHECK_INTERVAL

# Infinite loop to keep checking the job status
while true; do
    # Initialize flag to check if all jobs have completed
    allCompleted=true

    # Get the list of job IDs for jobs with names containing "GCM", "spaCRT", or "score"
    jobIds=$(qstat | grep -E "$searchPatterns" | awk '{print $1}')
    
    for jobId in $jobIds; do
    # Use qstat to check if the job is still listed
    jobStatus=$(qstat | grep "$jobId")

    if [ -n "$jobStatus" ]; then
        echo "Job $jobId containing one of the specified strings is still running."
        allCompleted=false
    else
        echo "Job $jobId containing one of the specified strings has finished."
    fi
    done

    # If all jobs have completed, exit the loop
    if $allCompleted; then
        echo "All jobs have completed."
        break
    else
        echo "Some jobs are still running. Checking again in $CHECK_INTERVAL seconds."
        sleep $CHECK_INTERVAL
    fi
done