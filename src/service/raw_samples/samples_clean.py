from samples_client import ClearFrequencyService

CFS = ClearFrequencyService()

CFS.cleanup_shm()

print("[CFS] Clean up complete\n")