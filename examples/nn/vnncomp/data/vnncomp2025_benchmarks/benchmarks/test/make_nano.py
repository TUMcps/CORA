'''
Stanley Bak
Make a nano onnx model that multiplies by 0.5
'''

import onnx
from onnx import TensorProto
from skl2onnx.algebra.onnx_ops import OnnxMatMul, OnnxAdd, OnnxRelu

import onnxruntime as ort

import numpy as np

def predict_with_onnxruntime(model_def, *inputs):
    'run an onnx model'
    
    sess = ort.InferenceSession(model_def.SerializeToString())
    names = [i.name for i in sess.get_inputs()]

    inp = dict(zip(names, inputs))
    res = sess.run(None, inp)
    names = [o.name for o in sess.get_outputs()]
    print(names, res)

    return res[0]

def main():
    'make test_nano.onnx'

    output_filename = "test_nano.onnx"

    input_name = 'model_in'
    inp_shape = (1,)
    a_mat = np.array([[0.5]], dtype=float)

    relu_node = OnnxRelu(OnnxMatMul(input_name, a_mat, output_names=['hidden']), output_names=['model_out'])
    #add_node = OnnxAdd(matmul_node, c_mats, output_names=[old_input_name])

    i = onnx.helper.make_tensor_value_info('i', TensorProto.FLOAT, inp_shape)
    #initial_types = [('float_input', FloatTensorType([None, 2]))]
        
    #model = convert_sklearn(matmul_node, initial_types=initial_types, target_opset=13)
    model = relu_node.to_onnx({input_name: i}, target_opset=14)

    print(model)

    onnx.checker.check_model(model)

    onnx.save_model(model, output_filename)
    print(f"Saved converted model to: {output_filename}")

    # run test
    print(f"f(1) = {predict_with_onnxruntime(model, [1.0])}")
    print(f"f(-1) = {predict_with_onnxruntime(model, [-1.0])}")

    
if __name__ == "__main__":
    main()
