# DisparityMap

Generates a depth map image as well as an image for ever step along the way from 2 stereo images:
 - SIFT to gather keypoints (completed)
 - Low thresholded Harris corner detector to filter edges (completed)
 - Taylor expansion around each keypoint for better accuracy (*optional)
 - 128 dimension SIFT feature vectors to match keypoints
 - Computes the fundamental matrix
 - RANSAC to obtain the best fundamental matrix (*optional)
 - Rectifies stereo images using the fundamental matrix
 - Computes disparity map from matching key points (completed)
 - Computes depth map from disparity map (depth = basline length * focal length / disparirty)

### Disclaimer
 - Matrix class by Github user "akalicki"
